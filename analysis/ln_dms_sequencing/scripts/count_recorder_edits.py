"""Count recorder editing per barcode for one RNA sample, with matched denominators.

Run from the repository root with the identifier of an RNA sample in the dataset manifest. One
streaming pass validates the construct, masks every measured base below the quality floor as N, and
tabulates each annotated adenosine over exactly the reads where that site is callable. A
masked base leaves only its own site's denominator; the rest of the read still counts.

Structure validation is editing-aware by default, because requiring an exact constant region
rejects precisely the edited molecules the assay measures. The strict pass rate is counted
alongside, so the size of that selection is reported rather than assumed.

"""

import os
import sys
import numpy as np
import pandas as pd
from analysis_helpers import (BASES, BASE_CODES, COUNT_DIRECTORY, LINKAGE_DIRECTORY, MASK, MIN_BASE_QUALITY,
                              PHRED_OFFSET, constant_regions, count_columns, encode, feature_checks, layout,
                              manifest, measured_sites, read_files, rejection_reason, require,
                              synchronized_fastq_records, verify_sites)

sample_id = sys.argv[1]
os.makedirs(COUNT_DIRECTORY, exist_ok=True)
row = manifest(role="rna").set_index("sample_id").loc[sample_id]
features = layout(sample_id)
verify_sites(sample_id)
paths = read_files(sample_id)
checks = feature_checks(features, roles=("required",))
require(len(checks) > 0, f"{sample_id} declares no required constant region")
regions = constant_regions(checks)

sites = measured_sites()
sites = sites.loc[sites.feature.isin(features)].reset_index(drop=True)
require(len(sites) > 0, f"{sample_id} measures no annotated sites")
site_ids = list(sites.site_id)
site_offsets = np.array([features[feature].offset + index for feature, index in zip(sites.feature, sites.feature_index)])
# Sites are gathered one read at a time rather than one site at a time, so a read contributes all
# of its measured bases in two vector operations.
site_groups = []
for read in sorted(set(features[feature].read for feature in sites.feature)):
    rows_on_read = np.array([index for index, feature in enumerate(sites.feature)
                             if features[feature].read == read])
    site_groups.append((read, site_offsets[rows_on_read], rows_on_read))
mask_code = ord(MASK)
site_rows = np.arange(len(sites))

linkage = pd.read_csv(f"{LINKAGE_DIRECTORY}/{row.linkage_dataset}.csv")
linkage = linkage.loc[linkage.status.eq("assigned")].reset_index(drop=True)
require(linkage.barcode.is_unique, "The linkage table assigns a barcode twice")
index_of = {barcode: index for index, barcode in enumerate(linkage.barcode)}
barcodes = len(index_of)

calls = np.zeros((barcodes, len(site_ids), len(BASES)), dtype=np.int64)
totals = np.zeros(barcodes, dtype=np.int64)
needed = {features[feature].read for feature in sites.feature} | {check.read for check in checks}
needed |= {features["barcode"].read}
reads = sorted(needed)
require(needed <= set(paths), f"{sample_id} does not declare every read its layout needs")

barcode_spec = features["barcode"]
barcode_alphabet = [ord(base) for base in barcode_spec.alphabet]
losses = {name: 0 for name in ["short_read", "barcode_low_quality", "barcode_alphabet", "unlinked"]}
losses.update({f"structure_{check.name}": 0 for check in checks})
total = strict_pass = aware_pass = accepted = 0

for record in synchronized_fastq_records({name: paths[name] for name in reads}):
    total += 1
    encoded = {name: (encode(sequence), encode(qualities).astype(np.int16) - PHRED_OFFSET)
               for name, (sequence, qualities) in record.items()}

    # Mask each measured base independently below the quality floor.
    measured = np.full(len(site_ids), mask_code, dtype=np.uint8)
    for read, offsets, rows_on_read in site_groups:
        bases, quality = encoded[read]
        if offsets.max() >= len(bases):
            continue
        measured[rows_on_read] = np.where(quality[offsets] >= MIN_BASE_QUALITY, bases[offsets], mask_code)
    columns = BASE_CODES[measured]

    strict = True
    aware = True
    reason = ""
    for group in regions:
        result = group.evaluate(*encoded[group.read])
        if result is None:
            reason = "short_read"
            break
        strict &= result[0]
        aware &= result[1]
    if not reason and not aware:
        reason = rejection_reason(checks, encoded)
    if reason:
        losses[reason] += 1
        continue
    aware_pass += 1
    strict_pass += strict

    bases, quality = encoded[barcode_spec.read]
    region = bases[barcode_spec.offset:barcode_spec.offset + barcode_spec.length]
    if len(region) < barcode_spec.length:
        losses["short_read"] += 1
        continue
    if (quality[barcode_spec.offset:barcode_spec.offset + barcode_spec.length] < MIN_BASE_QUALITY).any():
        losses["barcode_low_quality"] += 1
        continue
    if not np.isin(region, barcode_alphabet).all():
        losses["barcode_alphabet"] += 1
        continue
    position = index_of.get(region.tobytes().decode())
    if position is None:
        losses["unlinked"] += 1
        continue
    accepted += 1

    totals[position] += 1
    calls[position, site_rows, columns] += 1


require(accepted > 0, f"No accepted reads for {sample_id}")
require(sum(losses.values()) + accepted == total, "Read accounting does not conserve records")

columns_out = {"barcode": linkage.barcode, "insert_num": linkage.insert_num, "reads": totals}
for position, site in enumerate(site_ids):
    for order, base in enumerate(BASES):
        columns_out[f"pos_{site}_{base}"] = calls[:, position, order]
table = pd.DataFrame(columns_out)
require(list(table) == ["barcode", "insert_num"] + count_columns(site_ids),
        "The count table does not match its declared columns")
for site in site_ids:
    require(table[[f"pos_{site}_{base}" for base in BASES]].sum(axis=1).eq(table.reads).all(),
            f"Per-base counts of {site} do not sum to accepted reads")
table.loc[table.reads.gt(0)].to_csv(f"{COUNT_DIRECTORY}/{sample_id}_reads.csv", index=False)

accounting = [dict(sample_id=sample_id, stage="records", reads=total, fraction=1.0),
              dict(sample_id=sample_id, stage="structure_pass_editing_aware", reads=aware_pass,
                   fraction=aware_pass / total),
              dict(sample_id=sample_id, stage="structure_pass_strict", reads=strict_pass,
                   fraction=strict_pass / total)]
for name, count in sorted(losses.items(), key=lambda item: -item[1]):
    accounting.append(dict(sample_id=sample_id, stage=f"lost_{name}", reads=count, fraction=count / total))
accounting.append(dict(sample_id=sample_id, stage="accepted_reads", reads=accepted, fraction=accepted / total))
pd.DataFrame(accounting).to_csv(f"{COUNT_DIRECTORY}/{sample_id}_accounting.csv", index=False)

print(f"{sample_id}: {total} records, structure editing-aware {aware_pass / total:.4f}, "
      f"strict {strict_pass / total:.4f}, accepted {accepted} ({accepted / total:.4f})")
for name, count in sorted(losses.items(), key=lambda item: -item[1]):
    if count:
        print(f"  lost {name}: {count / total:.4f}")
