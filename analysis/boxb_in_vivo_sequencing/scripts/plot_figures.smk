"""Plotting workflow for the in vivo boxB RNA libraries. Run from the repository root.

The six figure families below form the output. FASTQ processing stays in the separate count
workflow, run_analysis.smk. Every rule here declares its count tables, annotations, script, and helpers, and its
outputs are exactly the PDFs and PNGs the output manifest registers for that analysis, so a figure
cannot be produced without being registered and a registered figure cannot go unbuilt. The export rule
then writes the plotted values of each numbered manuscript panel to the repository-root source_data/.
"""

import glob

import pandas as pd

ANALYSIS = "analysis/boxb_in_vivo_sequencing"
SCRIPTS = f"{ANALYSIS}/scripts"
ANNOTATIONS = f"{ANALYSIS}/annotations"
DATA = f"{ANALYSIS}/data"
FIGURES = f"{ANALYSIS}/figures"
IN_VITRO = "analysis/boxb_in_vitro_sequencing"
CONTAINER = "docker://ghcr.io/rasilab/r_python:2.4.2"
PRODUCTION_STAGE = "quality"

samples = pd.read_csv(f"{ANNOTATIONS}/sample_annotations.csv", comment="#")
SAMPLES = samples.sample_id.tolist()
comparisons = pd.read_csv(f"{ANNOTATIONS}/sample_comparisons.csv", comment="#")
outputs = pd.read_csv(f"{ANNOTATIONS}/figure_outputs.csv", comment="#")
assert not outputs.source.duplicated().any(), "The output manifest repeats a source"
assert not outputs.filename.duplicated().any(), "The output manifest repeats a filename"
assert outputs.filename.str.endswith(".pdf").all(), "Register every figure under its PDF filename"

COUNTS = expand(f"{DATA}/recorder_counts/{PRODUCTION_STAGE}/{{sample_id}}.csv.gz", sample_id=SAMPLES)
ANNOTATION_FILES = [f"{ANNOTATIONS}/{name}.csv" for name in
                    ["sample_annotations", "sample_info", "sample_comparisons", "barcode_annotations",
                     "recorder_positions", "recorder_backgrounds", "recorder_joint_events",
                     "boxb_structure_regions", "boxb_substitution_measurements", "boxb_wildtype_mutant_list",
                     "excluded_constructs", "figure_outputs"]]
HELPERS = [f"{SCRIPTS}/analysis_helpers.py", f"{SCRIPTS}/plot_style.py", f"{SCRIPTS}/in_vitro_reuse.py",
           f"{IN_VITRO}/scripts/analysis_helpers.py", f"{IN_VITRO}/scripts/plot_style.py"]
edit_count_panels = pd.read_csv(f"{ANNOTATIONS}/recorder_edit_count_panels.csv", comment="#")
edit_count_samples = comparisons.merge(
    edit_count_panels.assign(sample_id=edit_count_panels.sample_ids.str.split(";")).explode("sample_id")
    [["comparison_set", "sample_id"]].drop_duplicates(),
    on=["comparison_set", "sample_id"], validate="one_to_one").sample_id.unique()
site_panels = pd.read_csv(f"{ANNOTATIONS}/recorder_site_panels.csv", comment="#")
site_panel_samples = comparisons.merge(
    site_panels.assign(sample_id=site_panels.sample_ids.str.split(";")).explode("sample_id")
    [["comparison_set", "sample_id"]],
    on=["comparison_set", "sample_id"], validate="one_to_one")
single_mutant_correlations = pd.read_csv(f"{ANNOTATIONS}/boxb_substitution_correlations.csv", comment="#")
single_mutant_samples = comparisons.merge(
    single_mutant_correlations[["comparison_set", "sample_id"]].drop_duplicates(),
    on=["comparison_set", "sample_id"], validate="one_to_one").sample_id.unique()
effect_heatmaps = pd.read_csv(f"{ANNOTATIONS}/boxb_effect_heatmaps.csv", comment="#")
effect_heatmap_samples = comparisons.merge(effect_heatmaps[["comparison_set", "sample_id"]].drop_duplicates(),
                                           on=["comparison_set", "sample_id"], validate="one_to_one").sample_id.unique()
mutation_distance = pd.read_csv(f"{ANNOTATIONS}/boxb_mutation_distance.csv", comment="#")
mutation_distance_samples = comparisons.merge(mutation_distance[["comparison_set", "sample_id"]].drop_duplicates(),
                                              on=["comparison_set", "sample_id"], validate="one_to_one").sample_id.unique()
read_coverage = pd.read_csv(f"{ANNOTATIONS}/boxb_read_coverage.csv", comment="#")
read_coverage_samples = comparisons.merge(read_coverage[["comparison_set", "sample_id"]].drop_duplicates(),
                                         on=["comparison_set", "sample_id"], validate="one_to_one").sample_id.unique()
in_vitro_comparison = pd.read_csv(f"{ANNOTATIONS}/boxb_in_vitro_comparison.csv", comment="#")
in_vitro_comparison_samples = comparisons.merge(
    in_vitro_comparison[["comparison_set", "sample_id"]].drop_duplicates(),
    on=["comparison_set", "sample_id"], validate="one_to_one").sample_id.unique()
# The comparison reads the purified-enzyme reactions through the in vitro helpers, so both reactions of
# the named set are declared here and the in vitro annotations are taken whole: that analysis decides
# which of them its effect calculation reads, and this workflow should not second-guess it.
in_vitro_sets = pd.read_csv(f"{IN_VITRO}/annotations/boxb_reaction_sets.csv", comment="#")
in_vitro_reactions = sorted(in_vitro_sets.loc[
    in_vitro_sets.reaction_set.isin(in_vitro_comparison.in_vitro_reaction_set), "sample_id"].unique())
# The manuscript panels whose plotted values are exported, and the libraries and figure-data tables the
# export reads. A panel may only name a figure this analysis registers.
source_data_panels = pd.read_csv(f"{ANNOTATIONS}/source_data_panels.csv", comment="#")
assert not source_data_panels.panel.duplicated().any(), "The source-data annotation repeats a panel"
assert len(source_data_panels.merge(outputs[["source", "filename"]], on=["source", "filename"])) \
    == len(source_data_panels), "A source-data panel names an unregistered figure"
source_data_samples = sorted(set(edit_count_samples) | set(single_mutant_samples)
                             | set(site_panel_samples.sample_id))
SOURCE_DATA_TABLES = [f"{DATA}/figure_data/{source}.csv" for source in
                      [read_coverage.source.item(), mutation_distance.source.item(), effect_heatmaps.source.item()]]
SOURCE_DATA = expand("source_data/{panel}.csv", panel=source_data_panels.panel)


def registered(analysis):
    """Every file one analysis writes: its registered PDFs and the PNG of each."""
    selected = outputs.loc[outputs.analysis.eq(analysis), "filename"].tolist()
    assert selected, f"No registered output for {analysis}"
    return [f"{FIGURES}/{name}" for name in selected] + \
           [f"{FIGURES}/{name[:-4]}.png" for name in selected]


localrules: all


rule all:
    input:
        *[path for analysis in outputs.analysis.unique() for path in registered(analysis)],
        SOURCE_DATA,


rule export_source_data:
    """One tidy CSV of plotted values per manuscript panel, in the repository-root source_data/."""
    input:
        counts = expand(f"{DATA}/recorder_counts/{PRODUCTION_STAGE}/{{sample_id}}.csv.gz",
                        sample_id=source_data_samples),
        annotations = ANNOTATION_FILES, helpers = HELPERS,
        selections = [f"{ANNOTATIONS}/source_data_panels.csv", f"{ANNOTATIONS}/recorder_edit_count_panels.csv"],
        figure_data = SOURCE_DATA_TABLES,
        in_vitro_outputs = f"{IN_VITRO}/annotations/figure_outputs.csv",
        script = f"{SCRIPTS}/export_source_data.py",
    output: SOURCE_DATA
    container: CONTAINER
    shell: "python {input.script}"


rule plot_recorder_edit_counts:
    input:
        counts = expand(f"{DATA}/recorder_counts/{PRODUCTION_STAGE}/{{sample_id}}.csv.gz", sample_id=edit_count_samples),
        annotations = ANNOTATION_FILES,
        selection = f"{ANNOTATIONS}/recorder_edit_count_panels.csv", helpers = HELPERS,
        in_vitro_outputs = f"{IN_VITRO}/annotations/figure_outputs.csv",
        script = f"{SCRIPTS}/plot_recorder_edit_counts.py",
    output:
        figures = registered("plot_recorder_edit_counts"),
        statistics = expand(f"{DATA}/figure_data/{{source}}_statistics.csv", source=edit_count_panels.source),
    container: CONTAINER
    shell: "python {input.script}"


rule plot_recorder_site_placement:
    input:
        counts = lambda wildcards: expand(f"{DATA}/recorder_counts/{PRODUCTION_STAGE}/{{sample_id}}.csv.gz",
            sample_id=site_panel_samples.loc[site_panel_samples.comparison_set.eq(wildcards.comparison_set), "sample_id"]),
        annotations = ANNOTATION_FILES, helpers = HELPERS,
        selection = f"{ANNOTATIONS}/recorder_site_panels.csv",
        in_vitro_annotations = [f"{IN_VITRO}/annotations/{name}.csv" for name in
                                ["figure_outputs", "recorder_edit_summary_panels", "recorder_context_selection"]],
        script = f"{SCRIPTS}/plot_recorder_site_placement.py",
    output:
        pdf = f"{FIGURES}/recorder_sites_and_placements_tada_lambdaN_{{comparison_set}}.pdf",
        png = f"{FIGURES}/recorder_sites_and_placements_tada_lambdaN_{{comparison_set}}.png",
        statistics = f"{DATA}/figure_data/site_placement_{{comparison_set}}_statistics.csv",
    wildcard_constraints:
        comparison_set = "|".join(comparisons.comparison_set.unique()),
    container: CONTAINER
    shell: "python {input.script} {wildcards.comparison_set}"


rule plot_boxb_site_correlation:
    input:
        counts = expand(f"{DATA}/recorder_counts/{PRODUCTION_STAGE}/{{sample_id}}.csv.gz", sample_id=single_mutant_samples),
        annotations = ANNOTATION_FILES, helpers = HELPERS,
        selection = f"{ANNOTATIONS}/boxb_substitution_correlations.csv",
        in_vitro_outputs = f"{IN_VITRO}/annotations/figure_outputs.csv",
        script = f"{SCRIPTS}/plot_boxb_site_correlation.py",
    output: registered("plot_boxb_site_correlation")
    container: CONTAINER
    shell: "python {input.script}"


rule plot_boxb_effect_heatmaps:
    input:
        counts = expand(f"{DATA}/recorder_counts/{PRODUCTION_STAGE}/{{sample_id}}.csv.gz", sample_id=effect_heatmap_samples),
        annotations = ANNOTATION_FILES, helpers = HELPERS,
        selection = f"{ANNOTATIONS}/boxb_effect_heatmaps.csv",
        in_vitro_outputs = f"{IN_VITRO}/annotations/figure_outputs.csv",
        script = f"{SCRIPTS}/plot_boxb_effect_heatmaps.py",
    output:
        figures = registered("plot_boxb_effect_heatmaps"),
        tables = expand(f"{DATA}/figure_data/{{source}}.csv", source=effect_heatmaps.source),
    container: CONTAINER
    shell: "python {input.script}"


rule plot_boxb_mutation_distance:
    input:
        counts = expand(f"{DATA}/recorder_counts/{PRODUCTION_STAGE}/{{sample_id}}.csv.gz", sample_id=mutation_distance_samples),
        annotations = ANNOTATION_FILES, helpers = HELPERS,
        selection = f"{ANNOTATIONS}/boxb_mutation_distance.csv",
        in_vitro_outputs = f"{IN_VITRO}/annotations/figure_outputs.csv",
        script = f"{SCRIPTS}/plot_boxb_mutation_distance.py",
    output:
        figures = registered("plot_boxb_mutation_distance"),
        tables = expand(f"{DATA}/figure_data/{{source}}.csv", source=mutation_distance.source),
    container: CONTAINER
    shell: "python {input.script}"


rule plot_boxb_read_coverage:
    input:
        counts = expand(f"{DATA}/recorder_counts/{PRODUCTION_STAGE}/{{sample_id}}.csv.gz", sample_id=read_coverage_samples),
        annotations = ANNOTATION_FILES, helpers = HELPERS,
        selection = f"{ANNOTATIONS}/boxb_read_coverage.csv",
        in_vitro_outputs = f"{IN_VITRO}/annotations/figure_outputs.csv",
        script = f"{SCRIPTS}/plot_boxb_read_coverage.py",
    output:
        figures = registered("plot_boxb_read_coverage"),
        tables = expand(f"{DATA}/figure_data/{{source}}{{suffix}}.csv", source=read_coverage.source,
                        suffix=["", "_constructs"]),
    container: CONTAINER
    shell: "python {input.script}"


rule plot_in_vitro_comparison:
    input:
        counts = expand(f"{DATA}/recorder_counts/{PRODUCTION_STAGE}/{{sample_id}}.csv.gz",
                        sample_id=in_vitro_comparison_samples),
        annotations = ANNOTATION_FILES, helpers = HELPERS,
        selection = f"{ANNOTATIONS}/boxb_in_vitro_comparison.csv",
        in_vitro_annotations = sorted(glob.glob(f"{IN_VITRO}/annotations/*.csv")),
        in_vitro_counts = expand(f"{IN_VITRO}/data/summary_stats_combined/{{sample_id}}.csv.gz",
                                 sample_id=in_vitro_reactions),
        script = f"{SCRIPTS}/plot_in_vitro_comparison.py",
    output: registered("plot_in_vitro_comparison") + expand("source_data/{panel}.csv", panel=in_vitro_comparison.source_data_panel)
    threads: 1
    resources: mem_mb=2048
    container: CONTAINER
    shell: "python {input.script}"
