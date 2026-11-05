"""Count recorder editing, boxB editing, and QC for one in vivo boxB RNA library in a single pass.

The reporter and UMI reads are read together, with matching read identifiers and valid FASTQ records
required. The reverse primer handle is checked at its designed position rather than anywhere in the
read: the amplicon ends at the handle and this 79-base read stops before it, so a handle appearing
earlier marks a truncated molecule.

Every measurement here is tabulated over exactly the molecules its own denominator counts: a single
site over the molecules callable at that site, a joint event over those callable at all of its sites,
and the edit-count bins over those whose whole recorder is callable. Run from the repository root.

Column conventions. The recorder columns keep the read (antisense) orientation the in vitro analysis
uses, so `pos_<i>_C` is an edited molecule at the i-th measured site and `num_<n>_C` a molecule with n
edits across all eight; `num_0_N` counts the fully callable molecules. The boxB columns instead use
sense coordinates and sense bases, `boxb_<position>_<base>` for positions 1-19, because that is how
boxB is annotated everywhere else; an A→G edit there appears as `boxb_<position>_G`.
"""

import argparse
import gzip
import itertools as it
import re
import sys
from collections import Counter, defaultdict

import pandas as pd

BASES = ["A", "C", "G", "T", "N"]
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
MIN_BASE_QUALITY = 30
PHRED_OFFSET = 33
QUALITY_FLOOR = chr(MIN_BASE_QUALITY + PHRED_OFFSET)
HANDLE = "AGCACAACA"
# The directory the counts and QC are written under, and the `stage` column of both QC tables.
PRODUCTION_STAGE = "quality"
# The failure stages are mutually exclusive and ordered: a read is attributed to the first one it hits,
# so the stages and the accepted reads sum to the raw record count.
STAGES = ["barcode_unrecognized", "premature_handle", "missing_handle", "flank_mismatch", "insert_not_acgt",
          "umi_invalid", "duplicate_molecule"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_fastq(path):
    """Yield (identifier, sequence, quality) and fail on a malformed or truncated record."""
    handle = gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")
    with handle:
        while True:
            title = handle.readline()
            if not title:
                return
            sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
            require(title.startswith("@") and plus.startswith("+") and quality, f"Malformed FASTQ record in {path}")
            sequence, quality = sequence.rstrip("\n"), quality.rstrip("\n")
            require(len(sequence) == len(quality), f"Sequence and quality differ in length in {path}")
            yield title[1:].split()[0], sequence, quality


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_id")
    parser.add_argument("counts_output")
    parser.add_argument("attrition_output")
    parser.add_argument("umi_output")
    parser.add_argument("--analysis-directory", default="analysis/boxb_in_vivo_sequencing")
    parser.add_argument("--max-reads", type=int, default=0, help="stop after this many raw reads; 0 reads all")
    arguments = parser.parse_args()

    annotations = f"{arguments.analysis_directory}/annotations"
    samples = pd.read_csv(f"{annotations}/sample_annotations.csv", comment="#").set_index("sample_id", verify_integrity=True)
    require(arguments.sample_id in samples.index, f"Unannotated library {arguments.sample_id}")
    sample = samples.loc[arguments.sample_id]
    barcodes = pd.read_csv(f"{annotations}/barcode_annotations.csv", comment="#")
    require(not barcodes.reverse_complement.duplicated().any(), "Duplicate construct barcodes")
    positions = pd.read_csv(f"{annotations}/recorder_positions.csv", comment="#")
    require(not positions.position_id.duplicated().any() and not positions.read_index.duplicated().any(), "Invalid recorder positions")
    joint_events = pd.read_csv(f"{annotations}/recorder_joint_events.csv", comment="#")
    backgrounds = pd.read_csv(f"{annotations}/recorder_backgrounds.csv", comment="#").set_index("g_depleted", verify_integrity=True)

    # The measured recorder offsets, ordered as the pos_<i> columns number them.
    site_order = positions.sort_values("read_index")
    site_offsets = site_order.read_index.tolist()
    require(site_order.editing_column.tolist() == [f"pos_{i}_C" for i in range(len(site_offsets))],
            "Recorder positions must number the pos_<i> columns in read order")
    site_index = {row.position_id: index for index, row in enumerate(site_order.itertuples())}
    joint_members = {}
    for event in joint_events.itertuples():
        members = event.sites.split(";")
        require(len(set(members)) == len(members) and set(members) <= set(site_index), f"Invalid joint sites for {event.position_id}")
        joint_members[event.position_id] = [site_index[member] for member in members]
    n_sites = len(site_offsets)
    boxb_length = int(barcodes.boxb_length.iloc[0])

    # Fixed per-construct geometry, read once so the inner loop only slices.
    offset = int(sample.barcode_offset)
    constructs = {}
    for row in barcodes.itertuples():
        recorder = backgrounds.loc[row.g_depleted, "recorder_read_sequence"]
        target_start = offset + int(row.target_start)
        boxb_start = offset + int(row.boxb_start)
        variable_start = offset + int(row.variable_start)
        handle_start = offset + int(row.handle_start)
        # Masked boxB positions are the randomized ones; everything else is a constant site whose
        # editing is measurable, which is what makes boxB self-editing visible.
        window = set(range(int(row.variable_start), int(row.variable_start) + int(row.variable_length)))
        constructs[row.reverse_complement] = dict(
            oligo_name=row.oligo_name, recorder_length=len(recorder),
            target_start=target_start, boxb_start=boxb_start, handle_start=handle_start,
            variable_start=variable_start, variable_length=int(row.variable_length),
            site_positions=[target_start + site for site in site_offsets],
            upstream=re.compile(row.upstream_seq), upstream_start=target_start - len(row.upstream_seq),
            downstream=re.compile(row.downstream_seq), downstream_start=target_start + len(recorder),
            downstream_length=len(row.downstream_seq),
            boxb_sense_positions=[(boxb_start + k, boxb_length - k) for k in range(boxb_length)
                                  if int(row.boxb_start) + k not in window],
        )
    require(len({construct["handle_start"] for construct in constructs.values()}) == 1, "Constructs disagree on the handle position")
    handle_start = next(iter(constructs.values()))["handle_start"]
    # The handle is inside a read only when the read is long enough to reach past it, which the
    # 79-base layout is not; a handle before its designed position marks a truncated molecule.
    handle_in_read = int(sample.read1_length) >= handle_start + len(HANDLE)

    umi_start, umi_length = int(sample.umi_start), int(sample.umi_length)

    counts = defaultdict(Counter)
    seen = set()
    families = Counter()
    stage_counts = Counter()
    umi_sequences = Counter()
    unrecognized_barcodes = Counter()
    raw_reads = 0

    # Only the two streams deposited in SRA are read, so the counts reproduce from the public runs.
    streams = [read_fastq(sample.read1_fastq), read_fastq(sample.umi_fastq)]
    for records in it.zip_longest(*streams):
        require(all(record is not None for record in records), "Read streams ended at different records")
        raw_reads += 1
        identifiers = {record[0] for record in records}
        require(len(identifiers) == 1, f"Desynchronized read streams at record {raw_reads}: {identifiers}")
        read, umi_read = records
        sequence, quality = read[1], read[2]

        barcode = sequence[offset:offset + int(sample.barcode_length)]
        construct = constructs.get(barcode)
        if construct is None:
            stage_counts["barcode_unrecognized"] += 1
            unrecognized_barcodes[barcode] += 1
            continue
        # Position-aware architecture: the handle marks truncation before its designed position and
        # confirms a complete amplicon at it. Rejecting it anywhere discards every intact long read.
        if HANDLE in sequence[:handle_start]:
            stage_counts["premature_handle"] += 1
            continue
        if handle_in_read and sequence[handle_start:handle_start + len(HANDLE)] != HANDLE:
            stage_counts["missing_handle"] += 1
            continue
        if not (construct["upstream"].fullmatch(sequence[construct["upstream_start"]:construct["target_start"]])
                and construct["downstream"].fullmatch(sequence[construct["downstream_start"]:construct["downstream_start"] + construct["downstream_length"]])):
            stage_counts["flank_mismatch"] += 1
            continue
        insert = sequence[construct["variable_start"]:construct["variable_start"] + construct["variable_length"]]
        if len(insert) != construct["variable_length"] or not set(insert) <= {"A", "C", "G", "T"}:
            stage_counts["insert_not_acgt"] += 1
            continue
        umi = umi_read[1][umi_start:umi_start + umi_length]
        if len(umi) != umi_length or not set(umi) <= {"A", "C", "G", "T"}:
            stage_counts["umi_invalid"] += 1
            continue

        key = (barcode, insert)
        record = counts[key]
        record["read_counts"] += 1
        umi_sequences[umi] += 1
        family = (barcode, insert, umi)
        if family in seen:
            families[family] += 1
            record["family_duplicate_reads"] += 1
            stage_counts["duplicate_molecule"] += 1
            continue
        seen.add(family)
        record["umi_counts"] += 1

        # Recorder bases below the quality floor become N, so a molecule leaves only the sites it
        # cannot call and stays in every other site's numerator and denominator.
        recorder = [sequence[position] if quality[position] >= QUALITY_FLOOR else "N"
                    for position in construct["site_positions"]]
        for index, base in enumerate(recorder):
            record[f"pos_{index}_{base}"] += 1
        if "N" not in recorder:
            for base in BASES:
                record[f"num_{recorder.count(base)}_{base}"] += 1
        for event, members in joint_members.items():
            if all(recorder[member] != "N" for member in members):
                # The read is antisense, so a read C is an edited (G) sense base and a read T unedited.
                state = "".join("G" if recorder[member] == "C" else "A" for member in members)
                record[f"joint_{event}_callable"] += 1
                record[f"joint_{event}_state_{state}"] += 1

        # boxB constant sites, reported in sense coordinates: an A→G edit inside boxB is what makes a
        # genotype call ambiguous, so it is measured rather than assumed absent.
        for position, sense_position in construct["boxb_sense_positions"]:
            base = sequence[position] if quality[position] >= QUALITY_FLOOR else "N"
            record[f"boxb_{sense_position}_{COMPLEMENT[base]}"] += 1

        insert_quality = quality[construct["variable_start"]:construct["variable_start"] + construct["variable_length"]]
        record["insert_q30_molecules"] += int(min(insert_quality) >= QUALITY_FLOOR)
        umi_quality = umi_read[2][umi_start:umi_start + umi_length]
        record["umi_q30_molecules"] += int(min(umi_quality) >= QUALITY_FLOOR)
        if arguments.max_reads and raw_reads >= arguments.max_reads:
            break

    require(raw_reads > 0, f"No reads for {arguments.sample_id}")
    accepted = sum(record["umi_counts"] for record in counts.values())
    require(raw_reads == accepted + sum(stage_counts[stage] for stage in STAGES),
            f"Raw reads do not partition into accepted molecules and failure stages for {arguments.sample_id}")

    columns = (["read_counts", "umi_counts", "family_duplicate_reads", "insert_q30_molecules", "umi_q30_molecules"]
               + [f"pos_{index}_{base}" for index, base in it.product(range(n_sites), BASES)]
               + [f"num_{number}_{base}" for number, base in it.product(range(n_sites + 1), BASES)]
               + [f"joint_{event}_callable" for event in joint_members]
               + [f"joint_{event}_state_{''.join(state)}" for event, members in joint_members.items()
                  for state in it.product("AG", repeat=len(members))]
               + [f"boxb_{position}_{base}" for position, base in it.product(range(1, boxb_length + 1), BASES)])
    table = pd.DataFrame([dict(barcode=key[0], insert=key[1], **{column: record[column] for column in columns})
                          for key, record in sorted(counts.items())])
    table.insert(0, "sample_id", arguments.sample_id)
    table = table.merge(barcodes[["reverse_complement", "oligo_name"]].rename(columns={"reverse_complement": "barcode"}),
                        on="barcode", validate="many_to_one")
    for event, members in joint_members.items():
        states = ["".join(state) for state in it.product("AG", repeat=len(members))]
        require(table[[f"joint_{event}_state_{state}" for state in states]].sum(axis=1).eq(table[f"joint_{event}_callable"]).all(),
                f"Joint state table does not sum to its callable coverage for {event}")
        edited = f"joint_{event}_state_{'G' * len(members)}"
        for index, member in enumerate(members):
            marginal = table[f"pos_{member}_C"]
            require(table[edited].le(marginal).all(), f"Joint edits exceed marginal {index} for {event}")
    require(table[[f"pos_{index}_{base}" for base in BASES]].sum(axis=1).eq(table.umi_counts).all() if n_sites else True,
            "Site base counts do not sum to the molecule count")
    for index in range(n_sites):
        require(table[[f"pos_{index}_{base}" for base in BASES]].sum(axis=1).eq(table.umi_counts).all(),
                f"Site {index} base counts do not sum to the molecule count")
    require(table[[f"num_{number}_C" for number in range(n_sites + 1)]].sum(axis=1).eq(table.num_0_N).all(),
            "Edit-count bins do not sum to the fully callable molecules")
    table.to_csv(arguments.counts_output, index=False, compression="gzip" if arguments.counts_output.endswith(".gz") else None)

    family_sizes = Counter(families.values())
    attrition = pd.DataFrame([dict(sample_id=arguments.sample_id, measure="raw_reads", kind="total", value=raw_reads)]
        + [dict(sample_id=arguments.sample_id, measure=stage, kind="exclusive_failure", value=stage_counts[stage]) for stage in STAGES]
        + [dict(sample_id=arguments.sample_id, measure=measure, kind="diagnostic", value=value) for measure, value in [
            ("accepted_molecules", accepted),
            ("accepted_reads", int(table.read_counts.sum())),
            ("recognized_constructs", int(table.barcode.nunique())),
            ("observed_inserts", len(table)),
            ("fully_callable_molecules", int(table.num_0_N.sum())),
            ("insert_q30_molecules", int(table.insert_q30_molecules.sum())),
            ("umi_q30_molecules", int(table.umi_q30_molecules.sum())),
            ("multi_read_families", len(families)),
            ("largest_family_reads", 1 + max(families.values(), default=0)),
            ("unrecognized_barcode_sequences", len(unrecognized_barcodes)),
            ("handle_inside_read1", int(handle_in_read)),
        ]]
        + [dict(sample_id=arguments.sample_id, measure=f"families_with_{size + 1}_reads", kind="family_size", value=count)
           for size, count in sorted(family_sizes.items())])
    attrition.insert(1, "stage", PRODUCTION_STAGE)
    attrition.to_csv(arguments.attrition_output, index=False)

    # UMI occupancy per construct and insert, against the capacity this library's UMI length allows.
    capacity = 4 ** umi_length
    occupancy = table[["sample_id", "barcode", "oligo_name", "insert", "umi_counts", "read_counts", "family_duplicate_reads"]].copy()
    occupancy.insert(1, "stage", PRODUCTION_STAGE)
    occupancy["umi_length"] = umi_length
    occupancy["umi_capacity"] = capacity
    occupancy["umi_occupancy"] = occupancy.umi_counts / capacity
    occupancy["reads_per_molecule"] = occupancy.read_counts / occupancy.umi_counts
    top = umi_sequences.most_common(1)[0] if umi_sequences else ("", 0)
    occupancy["library_most_frequent_umi"] = top[0]
    occupancy["library_most_frequent_umi_reads"] = top[1]
    occupancy["library_distinct_umi_sequences"] = len(umi_sequences)
    occupancy.to_csv(arguments.umi_output, index=False)

    print(f"[{PRODUCTION_STAGE}] {arguments.sample_id}: {raw_reads:,} raw reads -> {accepted:,} molecules over "
          f"{table.barcode.nunique()} constructs and {len(table):,} inserts")
    print(f"  maximum construct/insert UMI occupancy: {occupancy.umi_occupancy.max():.3f} of {capacity:,}")
    print(f"  most frequent UMI sequence: {top[0]} in {top[1]:,} accepted reads")
    for stage in STAGES:
        print(f"  {stage:<24} {stage_counts[stage]:>12,} ({100 * stage_counts[stage] / raw_reads:5.2f}%)")


if __name__ == "__main__":
    sys.exit(main())
