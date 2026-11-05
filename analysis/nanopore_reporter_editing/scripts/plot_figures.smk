"""Draw manuscript Figure 1C from the reporter read counts."""

ANALYSIS = "analysis/nanopore_reporter_editing"
CONTAINER = "docker://ghcr.io/rasilab/r_python:2.4.2"

rule all:
  input:
    f"{ANALYSIS}/figures/reporter_editing.pdf",
    f"{ANALYSIS}/figures/reporter_editing.png",
    "source_data/figure_1c.csv",
    f"{ANALYSIS}/tables/figure_1c_statistics.csv",

rule plot_reporter_editing:
  input:
    counts = f"{ANALYSIS}/source_data/tada_lambdaN_replicates_tada.csv",
    annotations = [f"{ANALYSIS}/annotations/figure_1c_samples.csv",
                   f"{ANALYSIS}/annotations/edit_classes.csv",
                   f"{ANALYSIS}/annotations/figure_outputs.csv",
                   "analysis/boxb_in_vitro_sequencing/annotations/recorder_edit_summary_panels.csv"],
    script = f"{ANALYSIS}/scripts/plot_reporter_editing.py",
  output:
    f"{ANALYSIS}/figures/reporter_editing.pdf",
    f"{ANALYSIS}/figures/reporter_editing.png",
    f"{ANALYSIS}/figures/reporter_editing.svg",
    "source_data/figure_1c.csv",
    f"{ANALYSIS}/tables/figure_1c_statistics.csv",
  container: CONTAINER
  shell: "python {input.script:q}"
