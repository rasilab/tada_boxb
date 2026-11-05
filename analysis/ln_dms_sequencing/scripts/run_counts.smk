"""Build barcode linkage and per-site read counts for the samples used in the final figures."""

import pandas as pd

ANALYSIS = "analysis/ln_dms_sequencing"
ANNOTATIONS = f"{ANALYSIS}/annotations"
SCRIPTS = f"{ANALYSIS}/scripts"
DATA = f"{ANALYSIS}/data"
CONTAINER = "docker://ghcr.io/rasilab/r_python:2.4.2"
HELPERS = f"{SCRIPTS}/analysis_helpers.py"
SHARED_ANNOTATIONS = [f"{ANNOTATIONS}/{name}.csv" for name in
                      ["dataset_manifest", "read_files", "read_layouts", "measured_adenosines", "insert_annotations"]]

manifest = pd.read_csv(f"{ANNOTATIONS}/dataset_manifest.csv").set_index("sample_id")
read_files = pd.read_csv(f"{ANNOTATIONS}/read_files.csv", comment="#")
layouts = pd.read_csv(f"{ANNOTATIONS}/read_layouts.csv")
sites = pd.read_csv(f"{ANNOTATIONS}/measured_adenosines.csv")
RNA = list(manifest.loc[manifest.role.eq("rna") & manifest.use_for_figures.eq("yes")].index)
if not RNA:
    raise ValueError("No RNA sample selected for the final figures")


def fastqs(wildcards):
    """Only the reads used by the annotated layout are pipeline dependencies."""
    features = layouts.loc[layouts.sample_id.eq(wildcards.sample_id)]
    needed = features.loc[features.role.isin(["required", "diagnostic"])
                          | features.feature.isin(["barcode", "insert", *sites.feature]), "read"].unique()
    selected = read_files.loc[read_files.sample_id.eq(wildcards.sample_id)
                             & read_files.read.isin(needed)]
    if not len(needed) or set(selected.read) != set(needed):
        raise ValueError(f"Missing FASTQ declarations for {wildcards.sample_id}")
    return list(selected.path)


localrules: all

rule all:
  input: expand(f"{DATA}/recorder_counts/{{sample_id}}_reads.csv", sample_id=RNA)

rule insert_annotations:
  """Enumerate the designed inserts that linkage and counting are annotated against."""
  input:
    script = f"{SCRIPTS}/generate_insert_annotations.py",
  output:
    f"{ANNOTATIONS}/insert_annotations.csv",
  log: f"{DATA}/manifests/generate_insert_annotations.log"
  container: CONTAINER
  shell: "python {input.script} &> {log}"

rule linkage:
  input:
    reads = fastqs,
    annotations = SHARED_ANNOTATIONS,
    script = f"{SCRIPTS}/build_linkage.py",
    helpers = HELPERS,
  output:
    table = f"{DATA}/barcode_linkage/{{sample_id}}.csv",
    accounting = f"{DATA}/barcode_linkage/{{sample_id}}_accounting.csv",
  log: f"{DATA}/barcode_linkage/{{sample_id}}_linkage.log"
  container: CONTAINER
  shell: "python {input.script} {wildcards.sample_id} &> {log}"

rule counts:
  input:
    reads = fastqs,
    linkage = lambda wildcards: f"{DATA}/barcode_linkage/{manifest.loc[wildcards.sample_id, 'linkage_dataset']}.csv",
    annotations = SHARED_ANNOTATIONS,
    script = f"{SCRIPTS}/count_recorder_edits.py",
    helpers = HELPERS,
  output:
    counts = f"{DATA}/recorder_counts/{{sample_id}}_reads.csv",
    accounting = f"{DATA}/recorder_counts/{{sample_id}}_accounting.csv",
  log: f"{DATA}/recorder_counts/{{sample_id}}_counts.log"
  container: CONTAINER
  shell: "python {input.script} {wildcards.sample_id} &> {log}"
