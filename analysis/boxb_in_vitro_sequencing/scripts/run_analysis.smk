"""Count-processing workflow; run from the repository root."""

import os
import re
import pandas as pd

ANALYSIS_DIRECTORY = "analysis/boxb_in_vitro_sequencing"
ANNOTATIONS = f"{ANALYSIS_DIRECTORY}/annotations"
SCRIPTS = f"{ANALYSIS_DIRECTORY}/scripts"
DATA = f"{ANALYSIS_DIRECTORY}/data"
COUNT_CONTAINER = "docker://ghcr.io/rasilab/python:1.0.0"
COMBINE_CONTAINER = "docker://ghcr.io/rasilab/r_python:2.4.2"
BARCODE_LENGTH = 10

samples = pd.read_csv(f"{ANNOTATIONS}/sample_annotations.csv", dtype=str).set_index("sample_id")
barcodes = pd.read_csv(f"{ANNOTATIONS}/barcode_annotations.csv", dtype=str).set_index("reverse_complement")

rule all:
    input:
        expand(f"{DATA}/split_fastq/{{sample_id}}/{{barcode}}_R{{read}}.fastq", sample_id=samples.index, barcode=barcodes.index, read=[1, 2]),
        expand(f"{DATA}/summary_stats/{{sample_id}}/{{barcode}}.csv", sample_id=samples.index, barcode=barcodes.index),
        expand(f"{DATA}/summary_stats_combined/{{sample_id}}.csv.gz", sample_id=samples.index),
        expand(f"{DATA}/boxb_editing_stats/{{sample_id}}/{{barcode}}.csv", sample_id=samples.index, barcode=barcodes.index),
        expand(f"{DATA}/boxb_stats_combined/{{sample_id}}.csv.gz", sample_id=samples.index)


def get_reads(wildcards):
    filenames = sorted(f"{DATA}/fastq/{name}" for name in os.listdir(f"{DATA}/fastq") if re.search(f"{wildcards.sample_id}_", name) and re.search("_R[123]_", name) and name.endswith(".fastq"))
    if not filenames:
        raise FileNotFoundError(f"No input FASTQ files for {wildcards.sample_id}")
    if len(filenames) > 3:
        raise ValueError(f"Too many input FASTQ files for {wildcards.sample_id}: {filenames}")
    return filenames


rule split_fastq_by_barcode:
    input:
        reads=get_reads,
        annotations=f"{ANNOTATIONS}/barcode_annotations.csv",
        script=f"{SCRIPTS}/split_by_barcode.awk"
    output:
        [f"{DATA}/split_fastq/{{sample_id}}/{barcode}_R1.fastq" for barcode in barcodes.index],
        [f"{DATA}/split_fastq/{{sample_id}}/{barcode}_R2.fastq" for barcode in barcodes.index]
    params:
        barcode_length=BARCODE_LENGTH,
        output_directory=f"{DATA}/split_fastq/{{sample_id}}"
    shell:
        "mkdir -p {params.output_directory:q}; awk -v barcode_len={params.barcode_length} -v output_dir={params.output_directory:q} -v read2={input.reads[1]:q} -f {input.script:q} {input.annotations:q} {input.reads[0]:q}"

rule calculate_summary_stats:
    input:
        R1=f"{DATA}/split_fastq/{{sample_id}}/{{barcode}}_R1.fastq",
        R2=f"{DATA}/split_fastq/{{sample_id}}/{{barcode}}_R2.fastq",
        annotations=f"{ANNOTATIONS}/barcode_annotations.csv",
        script=f"{SCRIPTS}/calculate_summary_stats_for_each_barcode.py"
    output: f"{DATA}/summary_stats/{{sample_id}}/{{barcode}}.csv"
    container: COUNT_CONTAINER
    shell: "python {input.script:q} {input.annotations:q} {input.R1:q} {input.R2:q} {output:q} {wildcards.barcode:q}"

rule combine_summary_stats:
    input:
        counts=[f"{DATA}/summary_stats/{{sample_id}}/{barcode}.csv" for barcode in barcodes.index],
        script=f"{SCRIPTS}/combine_barcode_summary_stats.py"
    output: f"{DATA}/summary_stats_combined/{{sample_id}}.csv.gz"
    container: COMBINE_CONTAINER
    shell: "python {input.script:q} {output:q} {wildcards.sample_id:q} {input.counts:q}"

rule calculate_boxb_loop_stats:
    input:
        R1=f"{DATA}/split_fastq/{{sample_id}}/{{barcode}}_R1.fastq",
        R2=f"{DATA}/split_fastq/{{sample_id}}/{{barcode}}_R2.fastq",
        annotations=f"{ANNOTATIONS}/barcode_annotations.csv",
        script=f"{SCRIPTS}/calculate_boxb_loop_editing.py"
    output: f"{DATA}/boxb_editing_stats/{{sample_id}}/{{barcode}}.csv"
    container: COUNT_CONTAINER
    shell: "python {input.script:q} {input.annotations:q} {input.R1:q} {input.R2:q} {output:q} {wildcards.barcode:q}"

rule combine_boxb_stats:
    input:
        counts=[f"{DATA}/boxb_editing_stats/{{sample_id}}/{barcode}.csv" for barcode in barcodes.index],
        script=f"{SCRIPTS}/combine_barcode_summary_stats.py"
    output: f"{DATA}/boxb_stats_combined/{{sample_id}}.csv.gz"
    container: COMBINE_CONTAINER
    shell: "python {input.script:q} {output:q} {wildcards.sample_id:q} {input.counts:q}"
