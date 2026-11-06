"""Derive every in vivo boxB construct coordinate from the design oligos and write the annotation.

The barcode offset, randomized window, recorder cassette, boxB region, placement distance and flank
patterns are all located inside each design oligo. The result is checked against the design, against
the in vitro annotation of the same constructs, and against the hairpin origins the loop counter
fixes. Run from the repository root.
"""

import re
import sys

import pandas as pd

ANALYSIS = "analysis/boxb_in_vivo_sequencing"
IN_VITRO = "analysis/boxb_in_vitro_sequencing"
DESIGN = f"{ANALYSIS}/annotations/oligo_design.csv"
BACKGROUNDS = f"{ANALYSIS}/annotations/recorder_backgrounds.csv"
REGIONS = f"{ANALYSIS}/annotations/boxb_structure_regions.csv"
IN_VITRO_BARCODES = f"{IN_VITRO}/annotations/barcode_annotations.csv"
OUTPUT = f"{ANALYSIS}/annotations/barcode_annotations.csv"

# The reverse primer that ends the amplicon. A full-length read reaches its first bases only when the
# read is longer than the amplicon; inside the amplicon it marks a prematurely truncated molecule.
HANDLE = "AGCACAACA"
FLANK_LENGTH = 10
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
EXPECTED_OLIGOS = 32
LAYOUT_PATTERN = r"^(spacer([35])_(\d+))"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def boxb_read_sequence():
    """The wild-type boxB as it is read: the construct is sequenced on the antisense strand."""
    regions = pd.read_csv(REGIONS).sort_values("start")
    rna = "".join(regions.rna_sequence)
    require(regions.iloc[0].start == 1 and regions.iloc[-1].end == len(rna), "Invalid boxB structure annotation")
    return reverse_complement(rna.replace("U", "T"))


def locate_boxb(read, variable_start, variable_length, boxb):
    """The one boxB placement that agrees with the design everywhere outside the randomized window."""
    window = range(variable_start, variable_start + variable_length)
    starts = [
        start
        for start in range(len(read) - len(boxb) + 1)
        if set(window) <= set(range(start, start + len(boxb)))
        and all(read[start + i] == base for i, base in enumerate(boxb) if start + i not in window)
    ]
    require(len(starts) == 1, "The boxB region is not uniquely placed by the design")
    return starts[0]


def flank_pattern(read, start, end, variable_start, variable_length):
    """The flank pattern: the designed bases, with a wildcard for each randomized base."""
    window = range(variable_start, variable_start + variable_length)
    return "".join("." if position in window else read[position] for position in range(start, end))


def derive():
    design = pd.read_csv(DESIGN)
    backgrounds = pd.read_csv(BACKGROUNDS, comment="#").set_index("g_depleted", verify_integrity=True)
    boxb = boxb_read_sequence()
    require(len(design) == EXPECTED_OLIGOS and not design.barcode.duplicated().any(), "Invalid oligo design")

    rows = []
    for oligo in design.itertuples():
        layout, side, spacer = re.match(LAYOUT_PATTERN, oligo.oligo_name).groups()
        barcode = reverse_complement(oligo.barcode)
        require(oligo.oligo_sequence.count(barcode) == 1, f"Barcode is not unique in {oligo.oligo_name}")
        barcode_offset = oligo.oligo_sequence.index(barcode)
        read = oligo.oligo_sequence[barcode_offset:]

        randomized = list(re.finditer(r"N+", read))
        require(len(randomized) == 1, f"Expected one randomized window in {oligo.oligo_name}")
        variable_start, variable_length = randomized[0].start(), len(randomized[0].group())

        g_depleted = "yes" if oligo.oligo_name.endswith("g_depleted") else "no"
        recorder = backgrounds.loc[g_depleted, "recorder_read_sequence"]
        require(read.count(recorder) == 1, f"Recorder is not unique in {oligo.oligo_name}")
        target_start = read.index(recorder)
        target_length = len(recorder)

        boxb_start = locate_boxb(read, variable_start, variable_length, boxb)
        window_start, window_end = (int(value) for value in oligo.oligo_name.split("boxb_random_")[1].split("_")[:2])
        require(variable_length == window_end - window_start + 1, f"Window width disagrees with its name in {oligo.oligo_name}")
        # The read is antisense, so a window ending at RNA position p starts at boxb_start + (19 - p).
        require(variable_start == boxb_start + len(boxb) - window_end, f"Window offset disagrees with the boxB origin in {oligo.oligo_name}")

        if target_start > boxb_start:
            require(side == "5", f"Placement side disagrees with the read order in {oligo.oligo_name}")
            spacing = target_start - (boxb_start + len(boxb))
            adjacent = read[boxb_start + len(boxb)]
        else:
            require(side == "3", f"Placement side disagrees with the read order in {oligo.oligo_name}")
            spacing = boxb_start - (target_start + target_length)
            adjacent = read[boxb_start - 1]
        require(adjacent == backgrounds.loc[g_depleted, "boxb_adjacent_base"], f"The base beside boxB is not the designed invariant in {oligo.oligo_name}")
        target_dist = spacing - 1
        # The spacer5 oligo names state the physical distance; the spacer3 names count how much spacer
        # was removed from the longest placement instead, so only the 5′ names can be checked against it.
        require(side != "5" or target_dist == int(spacer), f"Placement distance disagrees with the oligo name in {oligo.oligo_name}")

        require(read.count(HANDLE) == 1, f"Reverse primer is not unique in {oligo.oligo_name}")
        handle_start = read.index(HANDLE)

        rows.append(dict(
            barcode=oligo.barcode, reverse_complement=barcode, oligo_name=oligo.oligo_name,
            oligo_sequence=oligo.oligo_sequence, n_variants=oligo.n_variants,
            variable_start=variable_start, variable_length=variable_length,
            target_start=target_start, target_length=target_length,
            boxb_start=boxb_start, boxb_length=len(boxb), handle_start=handle_start,
            amplicon_length=len(read), barcode_offset_in_oligo=barcode_offset,
            layout=layout, window_start=window_start, window_end=window_end, designed_spacer_label=int(spacer),
            target_pos_to_boxb=int(side), target_dist=target_dist,
            variable_type="boxb", variable_subpos=f"{window_start}_{window_end}", g_depleted=g_depleted,
            upstream_seq=flank_pattern(read, target_start - FLANK_LENGTH, target_start, variable_start, variable_length),
            downstream_seq=flank_pattern(read, target_start + target_length, target_start + target_length + FLANK_LENGTH, variable_start, variable_length),
        ))

    table = pd.DataFrame(rows)
    require(table.barcode_offset_in_oligo.nunique() == 1 and table.amplicon_length.nunique() == 1, "Constructs disagree on the read origin or amplicon length")
    require(table.handle_start.eq(table.amplicon_length - len(HANDLE) - len("ACGAAGCCA")).all(), "The reverse primer does not end the amplicon")
    require(table.groupby("layout").boxb_start.nunique().eq(1).all(), "A layout carries two boxB origins")
    require(table.groupby("layout").target_dist.nunique().eq(1).all(), "A layout carries two placements")
    require(sorted(table.target_dist.unique()) == [0, 10, 20, 30], "Unexpected set of recorder placements")
    return table.sort_values(["variable_subpos", "g_depleted", "layout"]).reset_index(drop=True)


def check_against_in_vitro(table):
    """The same constructs are annotated for the in vitro library; the shared geometry must agree."""
    shared = ["variable_start", "variable_length", "target_start", "target_length", "target_pos_to_boxb", "target_dist", "variable_subpos", "g_depleted"]
    in_vitro = pd.read_csv(IN_VITRO_BARCODES)
    in_vitro = in_vitro.loc[in_vitro.variable_type.eq("boxb")]
    merged = table.merge(in_vitro, on="oligo_name", suffixes=("", "_in_vitro"), validate="one_to_one")
    require(len(merged) == len(table), "The in vitro annotation does not cover every in vivo construct")
    for column in shared:
        mismatched = merged.loc[merged[column].astype(str).ne(merged[f"{column}_in_vitro"].astype(str)), ["oligo_name", column, f"{column}_in_vitro"]]
        require(mismatched.empty, f"Derived {column} disagrees with the in vitro annotation:\n{mismatched.to_string(index=False)}")
    # The in vitro and in vivo pools reuse barcode strings for different constructs, so the barcode
    # itself must never be the join key between the two analyses.
    print(f"Barcode strings shared with the in vitro pool but naming another construct: "
          f"{merged.reverse_complement.ne(merged.reverse_complement_in_vitro).sum()} of {len(merged)}")


if __name__ == "__main__":
    derived = derive()
    check_against_in_vitro(derived)
    if "--dry-run" not in sys.argv:
        # No header comment: the columns are documented in annotations/README.md, so that GitHub
        # renders the file as a table.
        derived.to_csv(OUTPUT, index=False)
        print(f"\nWrote {len(derived)} constructs to {OUTPUT}")
