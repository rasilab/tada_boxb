#!/usr/bin/env python
"""
Calculate BoxB loop editing statistics from FASTQ sequencing data.

This script processes paired-end FASTQ files to analyze BoxB editing patterns
in specific loop regions, tracking UMI counts and nucleotide compositions.
"""

import pandas as pd
import re
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict
import itertools as it
import sys


def extract_insert_sequence(target, oligo_name, variable_start, variable_length, read_sequence):
    """Extract insert sequence based on oligo type."""
    target_5_var_positions = [8, 10, 11, 14, 16]
    target_3_var_positions = [0, 3, 5, 6, 8]
    
    if "boxb_random" in oligo_name:
        return read_sequence[variable_start:variable_start + variable_length]
    elif "target_random_5" in oligo_name:
        return ''.join(target[pos] for pos in target_5_var_positions)
    elif "target_random_3" in oligo_name:
        return ''.join(target[pos] for pos in target_3_var_positions)
    else:
        raise ValueError(f"Unknown oligo type: {oligo_name}")


def initialize_count_dict(loop_length, measured_length):
    """Initialize count dictionary with all possible keys."""
    count_dict = {'umi_counts': 0}
    nucleotides = ['A', 'G', 'T', 'C', 'N']
    
    # Initialize nucleotide frequency counts, which the edit-count bins take over the loop alone
    for num, nt in it.product(range(loop_length + 1), nucleotides):
        count_dict[f'num_{num}_{nt}'] = 0
    
    # Initialize position-specific counts, one per measured position
    for pos, nt in it.product(range(measured_length), nucleotides):
        count_dict[f'pos_{pos}_{nt}'] = 0
    
    return count_dict


def main():
    # Parse command line arguments
    annotations_file = sys.argv[1]
    fastq_file = sys.argv[2]
    fastq_umi = sys.argv[3]
    output_file = sys.argv[4]
    barcode = sys.argv[5]
    
    # Define constants
    UMI_START = 0
    UMI_LENGTH = 7
    SPACER3_BOXB_START = 60
    SPACER5_BOXB_START = 12
    BOXB_LENGTH = 19
    BOXB_LOOP_POSITIONS = [7, 9, 10]
    # Per-position counts also cover the adenosine that closes the stem, sense position 13 at read
    # offset 6. It is appended so pos_0 to pos_2 keep the loop positions they have always named.
    BOXB_MEASURED_POSITIONS = BOXB_LOOP_POSITIONS + [6]
    TARGET_EDIT_POSITIONS = [1, 2, 4, 7, 9, 12, 13, 15]
    TARGET_SEQUENCE = 'ATTATGGTGTAATTCTA'
    PRIMER_SEQUENCE = "AGCACAACA"
    MIN_BASE_QUALITY = 30
    PHRED_OFFSET = 33
    NUCLEOTIDES = ['A', 'G', 'T', 'C', 'N']

    # Load barcode annotations
    annotations = pd.read_csv(annotations_file).set_index('reverse_complement')
    annotations_barcode = annotations.loc[barcode]
    
    # Extract sequence validation parameters
    proper_upstream = annotations_barcode['upstream_seq']
    proper_downstream = annotations_barcode['downstream_seq']
    if pd.isna(proper_downstream):
        proper_downstream = ""
    
    # Pre-calculate extraction bounds
    target_start = annotations_barcode['target_start']
    target_length = annotations_barcode['target_length']
    target_end = target_start + target_length
    upstream_start = target_start - 10
    downstream_start = target_start + 17
    downstream_end = downstream_start + 10
    
    # Determine BoxB extraction parameters
    is_spacer3 = annotations_barcode['oligo_name'].startswith("spacer3")
    boxb_start = SPACER3_BOXB_START if is_spacer3 else SPACER5_BOXB_START
    boxb_end = boxb_start + BOXB_LENGTH
    
    # Initialize data structures
    count_table = defaultdict(lambda: initialize_count_dict(len(BOXB_LOOP_POSITIONS), len(BOXB_MEASURED_POSITIONS)))
    umi_sets = defaultdict(set)
    
    # Process paired FASTQ files
    with open(fastq_file, 'r') as f1, open(fastq_umi, 'r') as f2:
        R1_file = FastqGeneralIterator(f1)
        R2_file = FastqGeneralIterator(f2)
        
        for read1, read2 in zip(R1_file, R2_file):
            read_sequence = read1[1]
            
            # Skip reads containing primer sequence
            if PRIMER_SEQUENCE in read_sequence:
                continue
            
            # Extract sequences for validation
            target = read_sequence[target_start:target_end]
            upstream_seq = read_sequence[upstream_start:target_start]
            downstream_seq = read_sequence[downstream_start:downstream_end]
            
            # Validate sequence structure
            upstream_match = re.fullmatch(proper_upstream, upstream_seq)
            downstream_match = (downstream_seq == '' and proper_downstream == '') or \
                             re.fullmatch(proper_downstream, downstream_seq)
            
            if not (upstream_match and downstream_match):
                continue
            
            # Extract UMI and insert sequence
            umi = read2[1][UMI_START:UMI_START + UMI_LENGTH]
            
            try:
                insert = extract_insert_sequence(
                    target, 
                    annotations_barcode['oligo_name'],
                    annotations_barcode['variable_start'],
                    annotations_barcode['variable_length'],
                    read_sequence
                )
            except ValueError as e:
                print(f"Error: {e}")
                sys.exit(1)
            
            # Skip if UMI already seen for this insert
            if umi in umi_sets[insert]:
                continue
            
            umi_sets[insert].add(umi)
            count_dict = count_table[insert]
            count_dict['umi_counts'] += 1
            
            # Extract the measured BoxB positions, masking calls below the quality floor as N
            boxb = read_sequence[boxb_start:boxb_end]
            boxb_quality = read1[2][boxb_start:boxb_end]
            boxb_called_bases = {
                pos: boxb[pos] if ord(boxb_quality[pos]) - PHRED_OFFSET >= MIN_BASE_QUALITY else 'N'
                for pos in BOXB_MEASURED_POSITIONS
            }

            # Edit-count bins span every loop position, so they take only fully callable molecules
            loop_sequence = ''.join(boxb_called_bases[pos] for pos in BOXB_LOOP_POSITIONS)
            if 'N' not in loop_sequence:
                for nt in NUCLEOTIDES:
                    count = loop_sequence.count(nt)
                    count_dict[f"num_{count}_{nt}"] += 1

            # Per-position counts keep every molecule; umi_counts - pos_<i>_N is each denominator
            for index, pos in enumerate(BOXB_MEASURED_POSITIONS):
                count_dict[f'pos_{index}_{boxb_called_bases[pos]}'] += 1
    
    # Convert to DataFrame and save
    count_table_df = pd.DataFrame.from_dict(dict(count_table), orient='index')
    count_table_df.to_csv(output_file, index_label='insert')


if __name__ == "__main__":
    main()