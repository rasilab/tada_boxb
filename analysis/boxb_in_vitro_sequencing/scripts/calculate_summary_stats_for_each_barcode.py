#!/usr/bin/env python
"""
Calculate summary statistics for nucleotide composition in target editing regions.

This script processes paired-end FASTQ files to analyze editing patterns in target sequences,
tracking UMI counts and nucleotide compositions across different barcode types.
"""

import pandas as pd
import re
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict
import itertools as it
import sys


def initialize_count_dict(sequence_length, joint_columns=()):
    """Initialize count dictionary with all possible keys."""
    count_dict = {'umi_counts': 0}
    count_dict.update({column: 0 for column in joint_columns})
    # A joint event needs only its own sites callable, so it carries its own denominator
    count_dict.update({f'{column}_callable': 0 for column in joint_columns})
    nucleotides = ['A', 'G', 'T', 'C', 'N']
    
    # Initialize nucleotide frequency counts
    for num, nt in it.product(range(sequence_length + 1), nucleotides):
        count_dict[f'num_{num}_{nt}'] = 0
    
    # Initialize position-specific counts
    for pos, nt in it.product(range(sequence_length), nucleotides):
        count_dict[f'pos_{pos}_{nt}'] = 0
    
    return count_dict


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


def main():
    # Parse command line arguments
    annotations_file = sys.argv[1]
    fastq_file = sys.argv[2]
    fastq_umi = sys.argv[3]
    output_file = sys.argv[4]
    barcode = sys.argv[5]
    measurements_file = sys.argv[6] if len(sys.argv) > 6 else None
    positions_file = sys.argv[7] if measurements_file else None
    
    # Define constants
    UMI_START = 0
    UMI_LENGTH = 7
    TARGET_EDIT_POSITIONS = [1, 2, 4, 7, 9, 12, 13, 15]
    TARGET_SEQUENCE = 'ATTATGGTGTAATTCTA'
    PRIMER_SEQUENCE = "AGCACAACA"
    MIN_BASE_QUALITY = 30
    PHRED_OFFSET = 33
    NUCLEOTIDES = ['A', 'G', 'T', 'C', 'N']
    VARIABLE_TARGET_LENGTH = len(TARGET_EDIT_POSITIONS)
    
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

    # Optional joint events use the same accepted, first-seen UMI as site counts.
    joint_events = {}
    if measurements_file:
        measurements = pd.read_csv(measurements_file)
        positions = pd.read_csv(positions_file).set_index('position_id', verify_integrity=True)
        for measurement in measurements.itertuples():
            sites = measurement.sites.split(';')
            if len(sites) == 1:
                continue
            if len(set(sites)) != len(sites):
                raise ValueError('Duplicate sites in joint editing selection')
            members = []
            for site in sites:
                _, index, edited_base = positions.loc[site, 'editing_column'].split('_')
                index = int(index)
                if TARGET_EDIT_POSITIONS[index] != target_length - int(site[1:]):
                    raise ValueError('Recorder site and read coordinates disagree')
                members.append((index, edited_base))
            joint_events[f'joint_{measurement.position_id}'] = members
    
    # Initialize data structures
    count_table = defaultdict(lambda: initialize_count_dict(VARIABLE_TARGET_LENGTH, joint_events))
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
            
            # Extract variable target positions, masking calls below the quality floor as N
            target_quality = read1[2][target_start:target_end]
            variable_target = ''.join(
                target[pos] if ord(target_quality[pos]) - PHRED_OFFSET >= MIN_BASE_QUALITY else 'N'
                for pos in TARGET_EDIT_POSITIONS
            )

            # A joint event is measured on the molecules where its own sites are callable
            for column, members in joint_events.items():
                if all(variable_target[pos] != 'N' for pos, _ in members):
                    count_dict[f'{column}_callable'] += 1
                    count_dict[column] += int(all(variable_target[pos] == base for pos, base in members))

            # Edit-count bins count edits across every site, so they need the whole recorder
            # callable. num_0_N is their denominator.
            if 'N' not in variable_target:
                for nt in NUCLEOTIDES:
                    count = variable_target.count(nt)
                    count_dict[f"num_{count}_{nt}"] += 1

            # Per-site counts keep every molecule; umi_counts - pos_<site>_N is each site's denominator
            for pos, char in enumerate(variable_target):
                count_dict[f'pos_{pos}_{char}'] += 1
    
    # Convert to DataFrame and save
    count_table_df = pd.DataFrame.from_dict(dict(count_table), orient='index')
    count_table_df.to_csv(output_file, index_label='insert')


if __name__ == "__main__":
    main()
