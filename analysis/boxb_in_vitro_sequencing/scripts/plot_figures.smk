"""Figure workflow with PNG companions for heat maps and compact editor panels; counts are external inputs."""

import pandas as pd

ANALYSIS_DIRECTORY = "analysis/boxb_in_vitro_sequencing"
ANNOTATIONS = f"{ANALYSIS_DIRECTORY}/annotations"
SCRIPTS = f"{ANALYSIS_DIRECTORY}/scripts"
FIGURES = f"{ANALYSIS_DIRECTORY}/figures"
COUNTS = f"{ANALYSIS_DIRECTORY}/data/summary_stats_combined"
CONTAINER = "docker://ghcr.io/rasilab/r_python:2.4.2"
OUTPUT_MANIFEST = f"{ANNOTATIONS}/figure_outputs.csv"
COMMON = [OUTPUT_MANIFEST, f"{SCRIPTS}/analysis_helpers.py", f"{SCRIPTS}/plot_style.py"]
REACTIONS = f"{ANNOTATIONS}/sample_info.csv"
BARCODES = f"{ANNOTATIONS}/barcode_annotations.csv"
POSITIONS = f"{ANNOTATIONS}/recorder_positions.csv"
PANELS = f"{ANNOTATIONS}/recorder_edit_summary_panels.csv"
SITES = f"{ANNOTATIONS}/recorder_edit_summary_sites.csv"
SELECTION = f"{ANNOTATIONS}/recorder_context_selection.csv"
STEMS = f"{ANNOTATIONS}/boxb_wildtype_mutant_list.csv"
COMPONENTS = f"{ANNOTATIONS}/recorder_reaction_components.csv"
CONTEXT_SELECTION = f"{ANNOTATIONS}/recorder_context_heatmap_selection.csv"
STRUCTURE = f"{ANNOTATIONS}/boxb_structure_regions.csv"
MEASUREMENTS = f"{ANNOTATIONS}/boxb_substitution_measurements.csv"
CORRELATIONS = f"{ANNOTATIONS}/boxb_substitution_correlations.csv"
WINDOW_SELECTION = f"{ANNOTATIONS}/boxb_multimutant_heatmap_selection.csv"
REACTION_SETS = f"{ANNOTATIONS}/boxb_reaction_sets.csv"
HAIRPIN_FIGURES = f"{ANNOTATIONS}/boxb_hairpin_site_figures.csv"
RECRUITMENT_COMPARISON = f"{ANNOTATIONS}/boxb_recruitment_comparison.csv"
PANEL_SELECTION = f"{ANNOTATIONS}/source_data_panels.csv"
TABLES = f"{ANALYSIS_DIRECTORY}/tables"
# Source data is per manuscript panel, so it sits beside the manuscript rather than inside one analysis.
SOURCE_DATA = "source_data"
# The `short_key` of the reaction carrying no TadA, which every per-adenosine figure opens on.
CONTROL_SAMPLE_KEY = "no_enzyme"

panels = pd.read_csv(PANELS)
samples = pd.read_csv(REACTIONS)
selection = pd.read_csv(SELECTION)
metadata = selection[["sample_id", "tada_type", "tada_conc", "condition", "panel_order"]].drop_duplicates().sort_values("panel_order")
comparison_ids = set()
for panel in panels.to_dict("records"):
    selected = samples.loc[samples.short_key.isin(panel["short_keys"].split(";")) & samples.condition.isin(panel["conditions"].split(";")) & samples.tada_conc.isin(panel["concentrations"].split(";"))]
    comparison_ids.update(selected.sample_id)
outputs = pd.read_csv(OUTPUT_MANIFEST)
# A repeated source or filename silently doubles a rule's outputs, so catch it here rather than
# letting Snakemake report it as a duplicate output pattern far from its cause.
assert not outputs.source.duplicated().any(), f"Duplicate figure sources: {sorted(outputs.source[outputs.source.duplicated()])}"
assert not outputs.filename.duplicated().any(), f"Duplicate figure filenames: {sorted(outputs.filename[outputs.filename.duplicated()])}"
concentrations = panels.iloc[0].concentrations.split(";")
distance_ids = samples.loc[samples.tada_type.isin(metadata.tada_type) & samples.tada_conc.isin(concentrations) & samples.condition.eq(metadata.iloc[0].condition), "sample_id"]
comparison_outputs = [f"{FIGURES}/{name}" for name in outputs.loc[outputs.analysis.eq("plot_enzyme_comparisons"), "filename"]]
distance_outputs = [f"{FIGURES}/{name}" for name in outputs.loc[outputs.analysis.eq("plot_recorder_site_distances"), "filename"]]
context_outputs = [f"{FIGURES}/{name}" for name in outputs.loc[outputs.analysis.eq("plot_recorder_contexts"), "filename"]]
# A figure declares its own PNG companion beside its layout, so a README embed does not depend on a
# geometry setting that happens to differ from one.
sites = pd.read_csv(SITES)
png_groups = set(panels.loc[panels.png_companion, "figure_group"]) | set(sites.loc[sites.png_companion, "figure_group"])
comparison_outputs += [f"{FIGURES}/{name[:-4]}.png" for name in outputs.loc[outputs.source.isin(png_groups), "filename"]]
distance_outputs += [name[:-4] + ".png" for name in distance_outputs]
context_outputs += [name[:-4] + ".png" for name in context_outputs]
mutant_outputs = [f"{FIGURES}/{name}" for name in outputs.loc[outputs.analysis.eq("plot_boxb_substitutions"), "filename"]]
window_outputs = [f"{FIGURES}/{name}" for name in outputs.loc[outputs.analysis.eq("plot_boxb_multimutant_heatmaps"), "filename"]]
window_outputs += [f"{FIGURES}/{name[:-4]}.png" for name in outputs.loc[outputs.analysis.eq("plot_boxb_multimutant_heatmaps"), "filename"]]
mutation_distance_outputs = [f"{FIGURES}/{name}" for name in outputs.loc[outputs.analysis.eq("plot_boxb_mutation_distance"), "filename"]]
# Both boxB families draw every concentration set, so each needs the counts of all of its reactions.
reaction_sets = pd.read_csv(REACTION_SETS)
boxb_ids = sorted(reaction_sets.sample_id.unique())
measurements = pd.read_csv(MEASUREMENTS)
heatmap_sources = {f"{row.heatmap_source}_{set_id}" for row in measurements.loc[measurements.heatmap_source.notna()].itertuples() for set_id in row.reaction_sets.split(";")}
mutant_outputs += [f"{FIGURES}/{name[:-4]}.png" for name in outputs.loc[outputs.source.isin(heatmap_sources), "filename"]]
mutant_outputs += [f"{FIGURES}/{name[:-4]}.png" for name in outputs.loc[outputs.source.eq("boxb_single_mutant_site_correlation_time"), "filename"]]
# The all-variant correlations are read on their own, so every one of them carries a PNG companion.
correlations = pd.read_csv(CORRELATIONS)
designed_sources = {f"{row.source}_{set_id}" for row in correlations.loc[correlations.variant_set.eq("designed")].itertuples() for set_id in row.reaction_sets.split(";")}
mutant_outputs += [f"{FIGURES}/{name[:-4]}.png" for name in outputs.loc[outputs.source.isin(designed_sources), "filename"]]
mutation_distance_outputs += [f"{FIGURES}/{name[:-4]}.png" for name in outputs.loc[outputs.source.eq("boxb_mutation_distance_editors"), "filename"]]
# The direct reaction read against each split-recruitment reaction, one panel per row of the selection.
recruitment_selection = pd.read_csv(RECRUITMENT_COMPARISON)
recruitment_ids = set(recruitment_selection.y_sample_id) | set(recruitment_selection.x_sample_id)
recruitment_ids = sorted(recruitment_ids)
recruitment_outputs = [f"{FIGURES}/{name}" for name in outputs.loc[outputs.analysis.eq("plot_boxb_recruitment_comparison"), "filename"]]
recruitment_outputs += [name[:-4] + ".png" for name in recruitment_outputs]
# One per-adenosine figure per editor, each reading the loop counts of the reactions its own panel
# annotation names, which include the no-TadA control every one of them opens on.
hairpin_ids = sorted(set(pd.read_csv(HAIRPIN_FIGURES).sample_id))
hairpin_site_outputs = [f"{FIGURES}/{name}" for name in outputs.loc[outputs.analysis.eq("plot_boxb_hairpin_sites"), "filename"]]
hairpin_site_outputs += [name[:-4] + ".png" for name in hairpin_site_outputs]
# Every panel of the manuscript this folder supplies gets one table of the values its figure draws.
source_data_outputs = [f"{SOURCE_DATA}/{panel}.csv" for panel in pd.read_csv(PANEL_SELECTION).panel]

rule all:
    input: comparison_outputs + distance_outputs + context_outputs + mutant_outputs + window_outputs + mutation_distance_outputs + recruitment_outputs + hairpin_site_outputs + source_data_outputs


rule plot_boxb_hairpin_sites:
    input:
        counts=[f"{ANALYSIS_DIRECTORY}/data/boxb_stats_combined/{sample}.csv.gz" for sample in hairpin_ids],
        annotations=[REACTIONS, BARCODES, STEMS, STRUCTURE, PANELS, HAIRPIN_FIGURES],
        helpers=COMMON,
        script=f"{SCRIPTS}/plot_boxb_hairpin_sites.py"
    output: hairpin_site_outputs
    threads: 1
    resources: mem_mb=1024
    container: CONTAINER
    shell: "python {input.script:q}"

rule plot_enzyme_comparisons:
    input:
        counts=[f"{COUNTS}/{sample}.csv.gz" for sample in sorted(comparison_ids)],
        annotations=[REACTIONS, BARCODES, POSITIONS, PANELS, SITES, STEMS, COMPONENTS],
        helpers=COMMON,
        script=f"{SCRIPTS}/plot_enzyme_comparisons.py"
    output: comparison_outputs
    threads: 1
    resources: mem_mb=2048
    container: CONTAINER
    shell: "python {input.script:q}"

rule plot_recorder_site_distances:
    input:
        counts=[f"{COUNTS}/{sample}.csv.gz" for sample in distance_ids],
        annotations=[REACTIONS, BARCODES, POSITIONS, PANELS, STEMS, SELECTION],
        helpers=COMMON,
        script=f"{SCRIPTS}/plot_recorder_site_distances.py"
    output:
        figures=distance_outputs,
        statistics=[f"{ANALYSIS_DIRECTORY}/data/figure_data/{source}_statistics.csv"
                    for source in outputs.loc[outputs.analysis.eq("plot_recorder_site_distances"), "source"]]
    threads: 1
    resources: mem_mb=2048
    container: CONTAINER
    shell: "python {input.script:q}"

rule plot_recorder_contexts:
    input:
        counts=[f"{COUNTS}/{sample}.csv.gz" for sample in pd.read_csv(CONTEXT_SELECTION).sample_id.unique()],
        annotations=[REACTIONS, BARCODES, POSITIONS, PANELS, SELECTION, CONTEXT_SELECTION],
        helpers=COMMON,
        script=f"{SCRIPTS}/plot_recorder_contexts.py"
    output: context_outputs
    threads: 1
    resources: mem_mb=2048
    container: CONTAINER
    shell: "python {input.script:q}"

rule plot_boxb_substitutions:
    input:
        counts=[f"{COUNTS}/{sample}.csv.gz" for sample in boxb_ids],
        annotations=[REACTIONS, BARCODES, POSITIONS, PANELS, REACTION_SETS, STRUCTURE, MEASUREMENTS, CORRELATIONS],
        helpers=COMMON,
        script=f"{SCRIPTS}/plot_boxb_substitutions.py"
    output:
        figures=mutant_outputs,
        statistics=f"{ANALYSIS_DIRECTORY}/data/figure_data/boxb_substitution_correlation_statistics.csv"
    threads: 1
    resources: mem_mb=2048
    container: CONTAINER
    shell: "python {input.script:q}"

rule plot_boxb_multimutant_heatmaps:
    input:
        counts=[f"{COUNTS}/{sample}.csv.gz" for sample in boxb_ids],
        annotations=[REACTIONS, BARCODES, POSITIONS, PANELS, REACTION_SETS, STRUCTURE, WINDOW_SELECTION],
        helpers=COMMON,
        script=f"{SCRIPTS}/plot_boxb_multimutant_heatmaps.py"
    output: window_outputs
    threads: 1
    resources: mem_mb=2048
    container: CONTAINER
    shell: "python {input.script:q}"

rule plot_boxb_mutation_distance:
    input:
        counts=[f"{COUNTS}/{sample}.csv.gz" for sample in boxb_ids],
        annotations=[REACTIONS, BARCODES, POSITIONS, PANELS, REACTION_SETS, STRUCTURE, WINDOW_SELECTION],
        helpers=COMMON,
        script=f"{SCRIPTS}/plot_boxb_mutation_distance.py"
    output: mutation_distance_outputs
    threads: 1
    resources: mem_mb=2048
    container: CONTAINER
    shell: "python {input.script:q}"

rule plot_boxb_recruitment_comparison:
    input:
        counts=[f"{COUNTS}/{sample}.csv.gz" for sample in recruitment_ids],
        annotations=[REACTIONS, BARCODES, POSITIONS, PANELS, REACTION_SETS, STRUCTURE, WINDOW_SELECTION, RECRUITMENT_COMPARISON],
        helpers=COMMON,
        script=f"{SCRIPTS}/plot_boxb_recruitment_comparison.py"
    output: recruitment_outputs
    threads: 1
    resources: mem_mb=2048
    container: CONTAINER
    shell: "python {input.script:q}"


rule export_source_data:
    input:
        counts=[f"{COUNTS}/{sample}.csv.gz" for sample in sorted(set(comparison_ids) | set(distance_ids) | set(boxb_ids)
                                                                 | set(pd.read_csv(CONTEXT_SELECTION).sample_id) | set(recruitment_ids))],
        loop_counts=[f"{ANALYSIS_DIRECTORY}/data/boxb_stats_combined/{sample}.csv.gz" for sample in hairpin_ids],
        annotations=[REACTIONS, BARCODES, POSITIONS, PANELS, SITES, SELECTION, STEMS, COMPONENTS, CONTEXT_SELECTION,
                     STRUCTURE, MEASUREMENTS, REACTION_SETS, WINDOW_SELECTION, PANEL_SELECTION, RECRUITMENT_COMPARISON, HAIRPIN_FIGURES],
        helpers=COMMON,
        script=f"{SCRIPTS}/export_source_data.py"
    output: source_data_outputs
    threads: 1
    resources: mem_mb=4096
    container: CONTAINER
    shell: "python {input.script:q}"
