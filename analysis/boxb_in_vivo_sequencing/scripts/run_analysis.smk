"""Count workflow for the three in vivo boxB RNA libraries. Run from the repository root.

Every FASTQ, annotation, script, and output is declared, and the samples come from the annotation
rather than from a directory glob. One single-pass counter reads each library's raw streams and
writes its counts and QC together, applying the position-aware handle check and masking recorder
bases below Phred 30.
"""

import pandas as pd

ANALYSIS = "analysis/boxb_in_vivo_sequencing"
SCRIPTS = f"{ANALYSIS}/scripts"
ANNOTATIONS = f"{ANALYSIS}/annotations"
DATA = f"{ANALYSIS}/data"
CONTAINER = "docker://ghcr.io/rasilab/r_python:2.4.2"
STAGE = "quality"

samples = pd.read_csv(f"{ANNOTATIONS}/sample_annotations.csv", comment="#").set_index("sample_id", verify_integrity=True)
SAMPLES = samples.index.tolist()

ANNOTATION_FILES = [f"{ANNOTATIONS}/{name}.csv" for name in
                    ["sample_annotations", "sample_info", "barcode_annotations", "recorder_positions",
                     "recorder_joint_events", "recorder_backgrounds"]]


def raw_streams(sample_id):
    """Every raw FASTQ the counter opens for one library, taken from the annotation."""
    sample = samples.loc[sample_id]
    roles = ["read1_fastq", "umi_fastq"]
    return [sample[role] for role in roles if not pd.isna(sample[role])]


# Every rule but the target aggregator runs in the container, so none of them is local.
localrules: all


rule all:
    input:
        f"{ANNOTATIONS}/barcode_annotations.csv",
        expand(f"{DATA}/recorder_counts/{STAGE}/{{sample_id}}.csv.gz", sample_id=SAMPLES),
        expand(f"{DATA}/qc/{STAGE}/{{sample_id}}_attrition.csv", sample_id=SAMPLES),
        expand(f"{DATA}/qc/{STAGE}/{{sample_id}}_umi_occupancy.csv", sample_id=SAMPLES),


rule derive_barcode_annotations:
    """Regenerate every construct coordinate from the design oligos."""
    input:
        design = f"{ANNOTATIONS}/oligo_design.csv",
        backgrounds = f"{ANNOTATIONS}/recorder_backgrounds.csv",
        regions = f"{ANNOTATIONS}/boxb_structure_regions.csv",
        in_vitro = "analysis/boxb_in_vitro_sequencing/annotations/barcode_annotations.csv",
        script = f"{SCRIPTS}/derive_barcode_annotations.py",
    output:
        f"{ANNOTATIONS}/barcode_annotations.csv"
    container: CONTAINER
    shell: "python {input.script}"


rule count_recorder_molecules:
    input:
        fastq = lambda wildcards: raw_streams(wildcards.sample_id),
        annotations = ANNOTATION_FILES,
        script = f"{SCRIPTS}/count_recorder_molecules.py",
    output:
        counts = f"{DATA}/recorder_counts/{STAGE}/{{sample_id}}.csv.gz",
        attrition = f"{DATA}/qc/{STAGE}/{{sample_id}}_attrition.csv",
        umi = f"{DATA}/qc/{STAGE}/{{sample_id}}_umi_occupancy.csv",
    container: CONTAINER
    shell:
        """
        python {input.script} {wildcards.sample_id} {output.counts} {output.attrition} {output.umi}
        """
