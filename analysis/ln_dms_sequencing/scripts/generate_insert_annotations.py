"""Enumerate the designed λN inserts that linkage and counting are annotated against.

Each of the 22 wild-type codons is replaced by each of the 64 codons in turn. The 22 wild-type
replacements produce the same DNA sequence, so 21 of them collapse and 1,387 distinct designs
remain, with the canonical wild type retained once at the first position. Row numbers follow that
deduplicated order and are the insert identifiers the rest of the pipeline uses.
Run from the repository root.
"""

import itertools
import pandas as pd
from Bio.Seq import Seq

# The λN segment of the recorder reporter, read in the sense orientation. Codons are enumerated in
# the T, C, A, G order of the standard genetic-code table, which fixes the row order of the output.
WILD_TYPE_INSERT = "ATGGACGCTCAGACTCGCCGCCGCGAGCGCCGTGCAGAAAAACAAGCCCAGTGGAAAGCCGCCAAC"
CODON_BASE_ORDER = "TCAG"
CODON_LENGTH = 3
EXPECTED_DESIGNS = 1387
OUTPUT = "analysis/ln_dms_sequencing/annotations/insert_annotations.csv"

codons = ["".join(bases) for bases in itertools.product(CODON_BASE_ORDER, repeat=CODON_LENGTH)]
wt_codons = [WILD_TYPE_INSERT[start:start + CODON_LENGTH]
             for start in range(0, len(WILD_TYPE_INSERT), CODON_LENGTH)]
if len(WILD_TYPE_INSERT) % CODON_LENGTH or len(codons) != 4 ** CODON_LENGTH:
    raise ValueError("The wild-type insert is not a whole number of codons")

designs = []
for codon_position, wt_codon in enumerate(wt_codons, start=1):
    start = (codon_position - 1) * CODON_LENGTH
    for mut_codon in codons:
        insert_seq = WILD_TYPE_INSERT[:start] + mut_codon + WILD_TYPE_INSERT[start + CODON_LENGTH:]
        designs.append(dict(wt_codon=wt_codon, codon_position=codon_position, mut_codon=mut_codon,
                            insert_seq=insert_seq, wt_aa=str(Seq(wt_codon).translate()),
                            # Translate by the standard code alone. Every codon here is internal:
                            # λN is a C-terminal fusion to TadA-XTEN, so even position 1 is read
                            # through rather than initiated at, and a leading TTG or CTG encodes
                            # leucine rather than acting as an initiator.
                            mut_aa=str(Seq(mut_codon).translate()),
                            insert_seq_rc=str(Seq(insert_seq).reverse_complement())))
designs = pd.DataFrame(designs).drop_duplicates("insert_seq").reset_index(drop=True)
designs["row_num"] = range(1, len(designs) + 1)
if len(designs) != EXPECTED_DESIGNS or not designs.insert_seq.is_unique:
    raise ValueError(f"Expected {EXPECTED_DESIGNS} distinct designed inserts, got {len(designs)}")
if (designs.insert_seq.str.len() != len(WILD_TYPE_INSERT)).any():
    raise ValueError("A designed insert changed length")

designs.to_csv(OUTPUT, index=False, lineterminator="\n")
print(f"Wrote {len(designs):,} designed inserts over {len(wt_codons)} codon positions to {OUTPUT}; "
      f"the canonical wild type is row {designs.loc[designs.insert_seq.eq(WILD_TYPE_INSERT), 'row_num'].item()}.")
