"""One tidy source-data CSV of plotted values per manuscript panel, written to source_data/.

Every value comes from the shared helpers the figure scripts call, or from the figure-data tables
those scripts already write, so nothing is recalculated here and no panel can disagree with its own
source data. Each percentage carries the two counts it was divided from, and the denominator is named
in its column, because the denominators differ: a single recorder site is measured over the molecules
callable at that site, `umi_counts` minus that site's N count, while an edit-count bin counts edits
across all eight sites and is measured over the molecules callable at every one of them, `num_0_N`.

The exported values are the plotted ones, not an interpretation of them. The genotype a molecule is
binned under is read from the same molecule whose recorder is measured, so the limits recorded in the
README and in the interpretability audit apply to these tables exactly as they apply to the panels.
Run from the repository root.
"""

import os

import numpy as np
import pandas as pd

from analysis_helpers import (EXPECTED_DESIGNED_SEQUENCES, EXPECTED_SINGLE_SUBSTITUTIONS, annotation,
                              boxb_wildtype, comparison_samples, edit_count_fractions,
                              matched_stem_pairs, recorder_positions, require, sequence_effects,
                              substitution_measurements)
from plot_style import palettes

OUTPUT_DIRECTORY = "source_data"
FIGURE_DATA = "analysis/boxb_in_vivo_sequencing/data/figure_data"
COMPARISON_SET = "original"
EDITOR_SAMPLE = "i79_p21"
EDIT_COUNT_PANELS = ["figure_5c", "figure_s8a"]
SITES = ["A8", "A2"]
CONSTRUCT_KEYS = ["sample_id", "barcode", "insert"]
EDIT_BIN_LABELS = {"edited_1": "1 edit", "edited_2plus": "2+ edits"}
CONSTRUCT_TYPES = {"wt": "wild type", "mut": "stem mutant"}
OBSERVATION_TYPES = {"sequence": "pooled designed sequence", "control": "wild-type construct"}
AGREEMENT_TOLERANCE = 1e-9

panels = annotation("source_data_panels")
registered = annotation("figure_outputs")
require(not panels.panel.duplicated().any() and not panels.source.duplicated().any(),
        "A source-data panel or figure source is repeated")
require(len(panels.merge(registered[["source", "filename"]], on=["source", "filename"])) == len(panels),
        "A source-data panel names a figure the output manifest does not register")
sources = panels.set_index("panel").source
os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)
labels, _, _ = palettes()
selected = comparison_samples(COMPARISON_SET)
editor = selected.loc[selected.sample_id.eq(EDITOR_SAMPLE)]
require(len(editor) == 1, "The main-figure panels draw one editor library")
enzyme = labels[editor.iloc[0].enzyme]
site_columns = {site: site.lower() for site in SITES}
mean_column = "mean_" + "_".join(site_columns.values()) + "_percent_edited"


def write_panel(panel, table):
    """Write one panel's plotted values, after checking it carries no missing measurement."""
    require(len(table) > 0 and table.notna().all().all(), f"Incomplete source data for {panel}")
    table.to_csv(f"{OUTPUT_DIRECTORY}/{panel}.csv", index=False)
    print(f"{panel} ({sources[panel]}): {len(table)} rows, {len(table.columns)} columns")


# Figure 5B: observed reads per designed boxB sequence, pooled over placements and matching windows.
coverage = pd.read_csv(f"{FIGURE_DATA}/{sources['figure_5b']}.csv")
require(coverage.sample_id.eq(EDITOR_SAMPLE).all() and len(coverage) == EXPECTED_DESIGNED_SEQUENCES,
        "The coverage figure data does not hold one library's designed sequences")
write_panel("figure_5b", pd.DataFrame({
    "library": coverage.sample_id,
    "enzyme": enzyme,
    "coverage_rank": coverage["rank"],
    "boxb_sequence_dna": coverage.hairpin,
    "substitutions_from_wild_type": coverage.n_mutations,
    "is_wild_type": coverage.n_mutations.eq(0),
    "reads": coverage.read_counts,
    "molecules": coverage.umi_counts,
    "contributing_constructs_designed": coverage.n_designed_constructs,
    "contributing_constructs_observed": coverage.n_observed_constructs,
    "recorder_placements": coverage.n_placements,
    "randomized_windows": coverage.n_windows}))

# Figures 5C and S8A: edit-count bars, one row per construct and bin behind each drawn bar. Both
# panels are rows of the same annotation and differ only in the libraries they name.
edit_count_panels = annotation("recorder_edit_count_panels").set_index("source", verify_integrity=True)
for panel in EDIT_COUNT_PANELS:
    definition = edit_count_panels.loc[sources[panel]]
    libraries = selected.set_index("sample_id").loc[definition.sample_ids.split(";")].reset_index()
    pairs = matched_stem_pairs(libraries, substitution_measurements(COMPARISON_SET).iloc[:1])
    bins = edit_count_fractions(libraries)
    constructs = pairs.drop_duplicates(CONSTRUCT_KEYS)[
        CONSTRUCT_KEYS + ["enzyme", "oligo_name", "variable_subpos", "layout", "target_pos_to_boxb",
                          "target_dist", "insert_type"]]
    counts = constructs.merge(bins[CONSTRUCT_KEYS + ["num_0_N", *EDIT_BIN_LABELS]],
                              on=CONSTRUCT_KEYS, validate="one_to_one")
    long = counts.melt(id_vars=[column for column in counts.columns if column not in EDIT_BIN_LABELS],
                       value_vars=list(EDIT_BIN_LABELS), var_name="bin_column", value_name="edited_molecules")
    long["edit_bin"] = long.bin_column.map(EDIT_BIN_LABELS)
    long["construct_type"] = long.insert_type.map(CONSTRUCT_TYPES)
    long["percent_edited"] = 100 * long.edited_molecules / long.num_0_N
    require(np.isfinite(long.percent_edited).all(), "A selected construct has no fully callable molecule")
    bars = long.groupby(["sample_id", "edit_bin", "construct_type"]).percent_edited.agg(
        bar_mean_percent="mean", bar_sem_percent="sem", bar_constructs="size").reset_index()
    long = long.merge(bars, on=["sample_id", "edit_bin", "construct_type"], validate="many_to_one")
    write_panel(panel, pd.DataFrame({
        "library": long.sample_id,
        "enzyme": long.enzyme.map(labels),
        "edit_bin": long.edit_bin,
        "construct_type": long.construct_type,
        "barcode": long.barcode,
        "oligo_name": long.oligo_name,
        "randomized_window": long.variable_subpos,
        "recorder_placement": long.layout,
        "recorder_side_of_boxb": long.target_pos_to_boxb,
        "recorder_distance_nt": long.target_dist,
        "boxb_window_insert": long["insert"],
        "edited_molecules": long.edited_molecules,
        "molecules_callable_at_all_eight_sites": long.num_0_N,
        "percent_edited": long.percent_edited,
        "bar_mean_percent": long.bar_mean_percent,
        "bar_sem_percent": long.bar_sem_percent,
        "bar_constructs": long.bar_constructs}).sort_values(
            ["library", "edit_bin", "construct_type", "recorder_placement", "randomized_window",
             "boxb_window_insert"]))

# Figure 5D: every recorder adenosine, one row per construct and site behind each drawn point.
positions = recorder_positions().sort_values("sense_position").reset_index(drop=True)
site_definitions = positions.assign(sites=positions.position_id, label=positions.position_id,
                                    context=positions.derived_context)[
    ["position_id", "sites", "label", "context", "editing_column"]]
matched = matched_stem_pairs(editor, site_definitions)
points = matched.drop_duplicates(["sample_id", "position_id", "barcode", "insert"]).copy()
require(points.callable_umis.gt(0).all(), "A construct and recorder site carries no callable molecule")
points["percent_edited"] = 100 * points.fraction_edited
points["construct_type"] = points.insert_type.map(CONSTRUCT_TYPES)
drawn = points.groupby(["position_id", "construct_type"]).percent_edited.agg(
    point_mean_percent="mean", point_sem_percent="sem", point_constructs="size").reset_index()
points = points.merge(drawn, on=["position_id", "construct_type"], validate="many_to_one").merge(
    positions[["position_id", "sense_position"]], on="position_id", validate="many_to_one")
write_panel("figure_5d", pd.DataFrame({
    "library": points.sample_id,
    "enzyme": enzyme,
    "recorder_site": points.position_id,
    "recorder_site_sense_position": points.sense_position,
    "recorder_site_context": points.context,
    "construct_type": points.construct_type,
    "barcode": points.barcode,
    "oligo_name": points.oligo_name,
    "randomized_window": points.variable_subpos,
    "recorder_placement": points.layout,
    "recorder_side_of_boxb": points.target_pos_to_boxb,
    "recorder_distance_nt": points.target_dist,
    "boxb_window_insert": points["insert"],
    "edited_molecules": points.edited_umis,
    "molecules_callable_at_site": points.callable_umis,
    "percent_edited": points.percent_edited,
    "point_mean_percent": points.point_mean_percent,
    "point_sem_percent": points.point_sem_percent,
    "point_constructs": points.point_constructs}).sort_values(
        ["recorder_site_sense_position", "construct_type", "recorder_placement", "randomized_window",
         "boxb_window_insert"], ascending=[False, True, True, True, True]))

# Figures 5E, 5F, and 5G all read the same pooled per-sequence site counts, calculated once.
measurements = substitution_measurements(COMPARISON_SET).set_index("position_id").loc[SITES].reset_index()
require(measurements.sites.eq(measurements.position_id).all(), "Average marginal sites, not joint events")
effects, sequence_constructs = sequence_effects(editor, measurements)
sequences = effects.loc[effects.position_id.eq(SITES[0]), ["hairpin", "n_mutations"]].reset_index(drop=True)
for site, prefix in site_columns.items():
    rows = effects.loc[effects.position_id.eq(site)].set_index("hairpin", verify_integrity=True)
    require(len(rows) == EXPECTED_DESIGNED_SEQUENCES, "A recorder site is missing a designed sequence")
    sequences[f"{prefix}_edited_molecules"] = sequences.hairpin.map(rows.pooled_edited_umis)
    sequences[f"{prefix}_molecules_callable_at_site"] = sequences.hairpin.map(rows.pooled_umis)
    sequences[f"{prefix}_percent_edited"] = 100 * sequences.hairpin.map(rows.pooled_fraction_edited)
require(sequences.notna().all().all(), "A designed sequence is missing a pooled site measurement")
sequences[mean_column] = sequences[[f"{prefix}_percent_edited" for prefix in site_columns.values()]].mean(axis=1)

# Figure 5G: the heat-map tiles and the editing-loss logo drawn beneath them.
_, wildtype_rna, _ = boxb_wildtype()
tiles = pd.read_csv(f"{FIGURE_DATA}/{sources['figure_5g']}.csv")
tiles["wild_type_base"] = tiles.boxb_position.map(lambda position: wildtype_rna[position - 1])
require(tiles.is_wildtype.eq(tiles.nucleotide.eq(tiles.wild_type_base)).all(),
        "The heat-map grid disagrees with the wild-type boxB sequence")
heatmap = tiles.merge(sequences, on="hairpin", validate="many_to_one")
require(len(heatmap) == len(tiles) and np.allclose(heatmap[mean_column], heatmap.mean_editing_percent,
                                                   atol=AGREEMENT_TOLERANCE),
        "The drawn heat-map tiles disagree with the pooled site counts")
write_panel("figure_5g", pd.DataFrame({
    "library": heatmap.sample_id,
    "enzyme": enzyme,
    "boxb_position": heatmap.boxb_position,
    "wild_type_base": heatmap.wild_type_base,
    "mutant_base": heatmap.nucleotide,
    "is_wild_type_base": heatmap.is_wildtype,
    "boxb_sequence_dna": heatmap.hairpin,
    "a8_edited_molecules": heatmap.a8_edited_molecules,
    "a8_molecules_callable_at_site": heatmap.a8_molecules_callable_at_site,
    "a8_percent_edited": heatmap.a8_percent_edited,
    "a2_edited_molecules": heatmap.a2_edited_molecules,
    "a2_molecules_callable_at_site": heatmap.a2_molecules_callable_at_site,
    "a2_percent_edited": heatmap.a2_percent_edited,
    "tile_mean_a8_a2_percent_edited": heatmap.mean_editing_percent,
    "position_mean_percent_edited": heatmap.position_mean_editing_percent,
    "logo_position_editing_loss_percent": heatmap.position_editing_change_percent,
    "logo_letter_share_of_position": heatmap.relative_editing,
    "logo_letter_height_percent": heatmap.logo_letter_change_percent}).sort_values(
        ["boxb_position", "mutant_base"]))

# Figure 5E: every designed non-WT sequence, including measured zero editing.
mutant_sites = sequences.loc[sequences.n_mutations.gt(0)].copy()
require(len(mutant_sites) == EXPECTED_DESIGNED_SEQUENCES - 1, "Incomplete all-variant correlation")
mutant_sites.insert(0, "enzyme", enzyme)
mutant_sites.insert(0, "library", EDITOR_SAMPLE)
write_panel("figure_5e", mutant_sites.rename(columns={"hairpin": "boxb_sequence_dna", "n_mutations": "substitutions_from_wild_type"}))

# Figure 5F: one plotted observation per violin point, pooled sequences above zero substitutions and
# the distinct wild-type constructs at zero.
violin = pd.read_csv(f"{FIGURE_DATA}/{sources['figure_5f']}.csv")
mutants = violin.loc[violin.measurement.eq("pooled_sequence") & violin.n_mutations.gt(0)].merge(
    sequences, on=["hairpin", "n_mutations"], validate="one_to_one")
require(len(mutants) == EXPECTED_DESIGNED_SEQUENCES - 1
        and np.allclose(mutants[mean_column], mutants.mean_editing_percent, atol=AGREEMENT_TOLERANCE),
        "The drawn mutant violins disagree with the pooled site counts")
controls = sequence_constructs.loc[sequence_constructs.n_mutations.eq(0) & sequence_constructs.included_in_ratio]
wildtypes = controls.drop_duplicates("barcode")[
    ["barcode", "insert", "oligo_name", "variable_subpos", "layout", "target_pos_to_boxb", "target_dist",
     "hairpin", "n_mutations"]].reset_index(drop=True)
for site, prefix in site_columns.items():
    rows = controls.loc[controls.position_id.eq(site)].set_index("barcode", verify_integrity=True)
    wildtypes[f"{prefix}_edited_molecules"] = wildtypes.barcode.map(rows.edited_umis)
    wildtypes[f"{prefix}_molecules_callable_at_site"] = wildtypes.barcode.map(rows.callable_umis)
    wildtypes[f"{prefix}_percent_edited"] = 100 * wildtypes.barcode.map(rows.edited_umis / rows.callable_umis)
wildtypes[mean_column] = wildtypes[[f"{prefix}_percent_edited" for prefix in site_columns.values()]].mean(axis=1)
wildtypes = wildtypes.merge(violin.loc[violin.measurement.eq("WT_construct"), ["barcode", "mean_editing_percent"]],
                            on="barcode", validate="one_to_one")
require(np.allclose(wildtypes[mean_column], wildtypes.mean_editing_percent, atol=AGREEMENT_TOLERANCE),
        "The drawn wild-type construct points disagree with their own site counts")
measurement_columns = [f"{prefix}_{quantity}" for prefix in site_columns.values()
                       for quantity in ["edited_molecules", "molecules_callable_at_site", "percent_edited"]]
# A pooled sequence has no single construct, so its construct columns are deliberately empty.
observations = pd.concat([
    pd.DataFrame({"observation": OBSERVATION_TYPES["sequence"],
                  "substitutions_from_wild_type": mutants.n_mutations,
                  "boxb_sequence_dna": mutants.hairpin, "barcode": "", "oligo_name": "",
                  "randomized_window": "", "recorder_placement": "", "recorder_side_of_boxb": "",
                  "recorder_distance_nt": "", "boxb_window_insert": "",
                  **{column: mutants[column] for column in [*measurement_columns, mean_column]}}),
    pd.DataFrame({"observation": OBSERVATION_TYPES["control"],
                  "substitutions_from_wild_type": wildtypes.n_mutations,
                  "boxb_sequence_dna": wildtypes.hairpin, "barcode": wildtypes.barcode,
                  "oligo_name": wildtypes.oligo_name, "randomized_window": wildtypes.variable_subpos,
                  "recorder_placement": wildtypes.layout,
                  "recorder_side_of_boxb": wildtypes.target_pos_to_boxb,
                  "recorder_distance_nt": wildtypes.target_dist,
                  "boxb_window_insert": wildtypes["insert"],
                  **{column: wildtypes[column] for column in [*measurement_columns, mean_column]}})],
    ignore_index=True)
groups = observations.groupby("substitutions_from_wild_type")[mean_column].agg(
    group_median_percent="median", group_observations="size").reset_index()
observations = observations.merge(groups, on="substitutions_from_wild_type", validate="many_to_one")
observations.insert(0, "enzyme", enzyme)
observations.insert(0, "library", EDITOR_SAMPLE)
write_panel("figure_5f", observations.sort_values(
    ["substitutions_from_wild_type", "boxb_sequence_dna", "barcode"]))

print(f"\nWrote {len(panels)} panel tables to {OUTPUT_DIRECTORY}/. Each percentage carries its own "
      "numerator and denominator; a single recorder site is measured over the molecules callable at "
      "that site and an edit-count bin over the molecules callable at all eight sites.")
