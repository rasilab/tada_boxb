"""Generate the six figures and the summary statistics from the count tables."""

import pandas as pd

ANALYSIS = "analysis/ln_dms_sequencing"
ANNOTATIONS = f"{ANALYSIS}/annotations"
SCRIPTS = f"{ANALYSIS}/scripts"
DATA = f"{ANALYSIS}/data"
FIGURES = f"{ANALYSIS}/figures"
TABLES = f"{ANALYSIS}/tables"
SOURCE_DATA = "source_data"
CONTAINER = "docker://ghcr.io/rasilab/r_python:2.4.2"
HELPERS = [f"{SCRIPTS}/analysis_helpers.py", f"{SCRIPTS}/plot_style.py"]
BIOLOGY = [f"{ANNOTATIONS}/{name}.csv" for name in
           ["dataset_manifest", "read_layouts", "measured_adenosines", "insert_annotations"]]
FIGURE_NAMES = [
    "variant_read_coverage",
    "amino_acid_site_correlation_A2_A15",
    "recorder_site_profile",
    "position_editing_profile",
    "variant_class_distribution",
    "amino_acid_editing_logo",
]
# The manuscript panel each figure is exported as.
PANEL_NAMES = ["figure_6b", "figure_6c", "figure_6d", "figure_6e", "figure_6f", "figure_6h"]

manifest = pd.read_csv(f"{ANNOTATIONS}/dataset_manifest.csv")
samples = manifest.loc[manifest.role.eq("rna") & manifest.use_for_figures.eq("yes")]
if len(samples) != 1:
    raise ValueError("Select exactly one RNA sample for the retained figures")
INPUT_COUNTS = f"{DATA}/recorder_counts/{samples.sample_id.item()}_reads.csv"
INPUT_LINKAGE = f"{DATA}/barcode_linkage/{samples.linkage_dataset.item()}.csv"

localrules: all

rule all:
  input:
    figures = expand(f"{FIGURES}/{{name}}.{{extension}}", name=FIGURE_NAMES, extension=["pdf", "png"]),
    source_data = expand(f"{SOURCE_DATA}/{{panel}}.csv", panel=PANEL_NAMES),
    statistics = f"{TABLES}/wildtype_statistics.csv",

rule compute_wildtype_statistics:
  input:
    counts = INPUT_COUNTS,
    linkage = INPUT_LINKAGE,
    annotations = BIOLOGY,
    script = f"{SCRIPTS}/compute_wildtype_statistics.py",
    helpers = f"{SCRIPTS}/analysis_helpers.py",
  output:
    f"{TABLES}/wildtype_statistics.csv",
  log: f"{DATA}/manifests/compute_wildtype_statistics.log"
  container: CONTAINER
  shell: "python {input.script} &> {log}"

rule export_source_data:
  input:
    counts = INPUT_COUNTS,
    linkage = INPUT_LINKAGE,
    annotations = BIOLOGY + [f"{ANNOTATIONS}/amino_acid_groups.csv"],
    script = f"{SCRIPTS}/export_source_data.py",
    helpers = HELPERS,
  output:
    expand(f"{SOURCE_DATA}/{{panel}}.csv", panel=PANEL_NAMES),
  log: f"{DATA}/manifests/export_source_data.log"
  container: CONTAINER
  shell: "python {input.script} &> {log}"


rule plot_variant_read_coverage:
  input:
    counts = INPUT_COUNTS,
    linkage = INPUT_LINKAGE,
    annotations = BIOLOGY,
    script = f"{SCRIPTS}/plot_variant_read_coverage.py",
    helpers = HELPERS,
  output:
    f"{FIGURES}/variant_read_coverage.pdf",
    f"{FIGURES}/variant_read_coverage.png",
  log: f"{DATA}/manifests/plot_variant_read_coverage.log"
  container: CONTAINER
  shell: "python {input.script} &> {log}"

rule plot_amino_acid_site_correlation:
  input:
    counts = INPUT_COUNTS,
    linkage = INPUT_LINKAGE,
    annotations = BIOLOGY,
    script = f"{SCRIPTS}/plot_amino_acid_site_correlation.py",
    helpers = HELPERS,
  output:
    f"{FIGURES}/amino_acid_site_correlation_A2_A15.pdf",
    f"{FIGURES}/amino_acid_site_correlation_A2_A15.png",
  log: f"{DATA}/manifests/plot_amino_acid_site_correlation.log"
  container: CONTAINER
  shell: "python {input.script} &> {log}"

rule plot_recorder_site_profile:
  input:
    counts = INPUT_COUNTS,
    linkage = INPUT_LINKAGE,
    annotations = BIOLOGY,
    script = f"{SCRIPTS}/plot_recorder_site_profile.py",
    helpers = HELPERS,
  output:
    f"{FIGURES}/recorder_site_profile.pdf",
    f"{FIGURES}/recorder_site_profile.png",
  log: f"{DATA}/manifests/plot_recorder_site_profile.log"
  container: CONTAINER
  shell: "python {input.script} &> {log}"

rule plot_position_editing_profile:
  input:
    counts = INPUT_COUNTS,
    linkage = INPUT_LINKAGE,
    annotations = BIOLOGY,
    script = f"{SCRIPTS}/plot_position_editing_profile.py",
    helpers = HELPERS,
  output:
    f"{FIGURES}/position_editing_profile.pdf",
    f"{FIGURES}/position_editing_profile.png",
  log: f"{DATA}/manifests/plot_position_editing_profile.log"
  container: CONTAINER
  shell: "python {input.script} &> {log}"

rule plot_variant_class_distribution:
  input:
    counts = INPUT_COUNTS,
    linkage = INPUT_LINKAGE,
    annotations = BIOLOGY,
    script = f"{SCRIPTS}/plot_variant_class_distribution.py",
    helpers = HELPERS,
  output:
    f"{FIGURES}/variant_class_distribution.pdf",
    f"{FIGURES}/variant_class_distribution.png",
  log: f"{DATA}/manifests/plot_variant_class_distribution.log"
  container: CONTAINER
  shell: "python {input.script} &> {log}"

rule plot_amino_acid_editing_logo:
  input:
    counts = INPUT_COUNTS,
    linkage = INPUT_LINKAGE,
    annotations = BIOLOGY + [f"{ANNOTATIONS}/amino_acid_groups.csv"],
    script = f"{SCRIPTS}/plot_amino_acid_editing_logo.py",
    helpers = HELPERS,
  output:
    f"{FIGURES}/amino_acid_editing_logo.pdf",
    f"{FIGURES}/amino_acid_editing_logo.png",
  log: f"{DATA}/manifests/plot_amino_acid_editing_logo.log"
  container: CONTAINER
  shell: "python {input.script} &> {log}"
