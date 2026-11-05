"""Rebuild the barcode-to-variant linkage of one library from its reads, with qualities.

Run from the repository root with the identifier of a linkage sample in the dataset manifest. Every
read is parsed with its base qualities; a base below the quality floor masks its own position, and a
barcode or insert carrying a masked base is rejected with its reason recorded rather than silently
dropped. Inserts are tabulated before any design filter, so a barcode competing with an off-design
sequence is visible instead of being hidden by restricting to the design first. Libraries are kept
separate: a barcode shared by two independently generated pools is not treated as a conflict.
"""

import os
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
from analysis_helpers import (EXCLUDE_NEIGHBOURS, LINKAGE_DIRECTORY, MAX_OFF_DESIGN_SHARE, MIN_BASE_QUALITY,
                              MIN_LINKAGE_PURITY, MIN_LINKAGE_SUPPORT, MIN_REPORTED_SUPPORT, PHRED_OFFSET,
                              constant_regions, designed_insert_keys, encode, encode_sequence, feature_checks,
                              hamming_one_groups, inserts, layout, manifest, packable, read_files,
                              rejection_reason, require, synchronized_fastq_records, verify_sites)


sample_id = sys.argv[1]
os.makedirs(LINKAGE_DIRECTORY, exist_ok=True)
row = manifest(role="linkage").set_index("sample_id").loc[sample_id]
features = layout(sample_id)
verify_sites(sample_id)
require("insert" in features and "barcode" in features, f"{sample_id} declares no insert or barcode")
checks = feature_checks(features, roles=("required",))
regions = constant_regions(checks)
diagnostics = feature_checks(features, roles=("diagnostic",))
designed = designed_insert_keys(features["insert"].orientation)
designs = inserts().set_index("row_num")

paths = read_files(sample_id)
needed = {features["barcode"].read, features["insert"].read} | {check.read for check in checks + diagnostics}
reads = sorted(needed)
require(needed <= set(paths), f"{sample_id} does not declare every read its layout needs")

barcode_spec = features["barcode"]
insert_spec = features["insert"]
barcode_alphabet = [ord(base) for base in barcode_spec.alphabet]
rejections = defaultdict(int)
diagnostic_hits = defaultdict(int)
support = defaultdict(lambda: defaultdict(int))
accepted = total = 0

for record in synchronized_fastq_records({name: paths[name] for name in reads}):
    total += 1
    encoded = {name: (encode(sequence), encode(qualities).astype(np.int16) - PHRED_OFFSET)
               for name, (sequence, qualities) in record.items()}
    reason = ""
    for group in regions:
        result = group.evaluate(*encoded[group.read])
        if result is None:
            reason = "short_read"
            break
        if not result[1]:
            reason = rejection_reason(checks, encoded)
            break
    if not reason:
        for check in diagnostics:
            result = check.evaluate(*encoded[check.read])
            if result is not None and result[1]:
                diagnostic_hits[check.name] += 1
    if not reason:
        bases, quality = encoded[barcode_spec.read]
        region = bases[barcode_spec.offset:barcode_spec.offset + barcode_spec.length]
        low = quality[barcode_spec.offset:barcode_spec.offset + barcode_spec.length] < MIN_BASE_QUALITY
        if len(region) < barcode_spec.length:
            reason = "short_read"
        elif low.any():
            reason = "barcode_low_quality"
        elif not np.isin(region, barcode_alphabet).all():
            reason = "barcode_alphabet"
        else:
            barcode = region.tobytes().decode()
    if not reason:
        bases, quality = encoded[insert_spec.read]
        region = bases[insert_spec.offset:insert_spec.offset + insert_spec.length]
        low = quality[insert_spec.offset:insert_spec.offset + insert_spec.length] < MIN_BASE_QUALITY
        if len(region) < insert_spec.length:
            reason = "short_read"
        elif low.any() or not packable(region):
            reason = "insert_low_quality"
        else:
            support[barcode][encode_sequence(region.tobytes().decode())] += 1
            accepted += 1
    if reason:
        rejections[reason] += 1


require(accepted > 0, f"No accepted linkage reads for {sample_id}")

# Collapse each barcode's inserts into its dominant assignment, keeping the competition visible.
#
# Purity measures competition between designed library members, which is what a barcode shared by
# two variants looks like. Off-design reads are counted and reported separately rather than folded
# into purity: a single sequencing error anywhere in the 66-nt insert creates a unique off-design
# sequence, so at this run's depth an any-difference purity would reject most barcodes for read
# error rather than for genuine ambiguity.
rows = []
for barcode, counts in support.items():
    total_support = sum(counts.values())
    designed_counts = {designed[key]: count for key, count in counts.items() if key in designed}
    designed_support = sum(designed_counts.values())
    if designed_support < MIN_REPORTED_SUPPORT:
        continue
    ranked = sorted(designed_counts.items(), key=lambda item: -item[1])
    dominant_insert, dominant_support = ranked[0]
    rows.append(dict(barcode=barcode, support=total_support, designed_support=designed_support,
                     off_design_support=total_support - designed_support,
                     distinct_inserts=len(counts), distinct_designed_inserts=len(designed_counts),
                     insert_num=dominant_insert, dominant_support=dominant_support,
                     runner_up_support=ranked[1][1] if len(ranked) > 1 else 0,
                     purity=dominant_support / designed_support,
                     off_design_share=(total_support - designed_support) / total_support))
table = pd.DataFrame(rows)
require(len(table) > 0, f"No barcode of {sample_id} reaches the reported support floor")

neighbours = hamming_one_groups(table.barcode)
table["has_neighbour"] = table.barcode.isin(neighbours)
table["status"] = "assigned"
table.loc[table.designed_support.lt(MIN_LINKAGE_SUPPORT), "status"] = "low_support"
table.loc[table.purity.lt(MIN_LINKAGE_PURITY) & table.status.eq("assigned"), "status"] = "ambiguous_insert"
table.loc[table.off_design_share.gt(MAX_OFF_DESIGN_SHARE) & table.status.eq("assigned"), "status"] = "off_design"
if EXCLUDE_NEIGHBOURS:
    table.loc[table.has_neighbour & table.status.eq("assigned"), "status"] = "barcode_neighbour"
table = table.merge(designs[["codon_position", "wt_aa", "mut_aa", "mut_codon", "variant_class", "is_wild_type"]],
                    left_on="insert_num", right_index=True, how="left")
table.sort_values(["support", "barcode"], ascending=[False, True]).to_csv(f"{LINKAGE_DIRECTORY}/{sample_id}.csv", index=False)

accounting = [dict(sample_id=sample_id, stage="records", reads=total, fraction=1.0)]
for reason, count in sorted(rejections.items(), key=lambda item: -item[1]):
    accounting.append(dict(sample_id=sample_id, stage=f"rejected_{reason}", reads=count, fraction=count / total))
for name, count in sorted(diagnostic_hits.items()):
    accounting.append(dict(sample_id=sample_id, stage=f"diagnostic_{name}", reads=count, fraction=count / total))
accounting.append(dict(sample_id=sample_id, stage="accepted", reads=accepted, fraction=accepted / total))
require(sum(rejections.values()) + accepted == total, "Linkage read accounting does not conserve records")
pd.DataFrame(accounting).to_csv(f"{LINKAGE_DIRECTORY}/{sample_id}_accounting.csv", index=False)

chosen = table.loc[table.status.eq("assigned")]
print(f"{sample_id}: {total} records, {accepted} accepted ({accepted / total:.3f})")
for reason, count in sorted(rejections.items(), key=lambda item: -item[1]):
    print(f"  rejected {reason}: {count / total:.4f}")
for name, count in sorted(diagnostic_hits.items()):
    print(f"  diagnostic {name}: {count / total:.4f}")
print(f"  reported barcodes {len(table)}, assigned {len(chosen)}, "
      f"designed variants {chosen.insert_num.nunique()}, "
      f"cells {len(chosen[['codon_position', 'mut_aa']].drop_duplicates().dropna())}")
print("  status counts:", table.status.value_counts().to_dict())
