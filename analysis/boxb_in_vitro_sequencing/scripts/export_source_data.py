"""One tidy CSV per manuscript panel, taken from the helpers that draw the figures themselves.

Nothing is recalculated here: every value comes from the shared calculations in `analysis_helpers`,
so a number in `source_data/` is the number its panel draws. Each file names its reaction, its
measurement, and the molecule counts the plotted percentage rests on. Which panel reads which figure
is declared in the panel annotation, which the plotting workflow reads for the same list.
"""

import itertools
import os
import numpy as np
import pandas as pd
from analysis_helpers import (CONTROL_ENZYME, annotation, comparison_data, context_data, designed_variant_effects, distance_data,
                              hairpin_site_data, hairpin_site_figures, SAMPLE_KEYS, sequence_effects, multimutant_window_data, mutation_distance_data, require,
                              selected_measurement, single_mutant_data)
from plot_style import palettes

SOURCE_DIRECTORY = "source_data"
MEAN_MEASUREMENT = "mean_A8_A2"
# The two recorder adenosines the boxB panels average, x axis first where a panel pairs them.
CORRELATION_SITES = ["A8", "A2"]
FAMILIES = ["edit_count", "site", "hairpin", "recruitment", "distance", "context", "heatmap", "designed", "substitution",
            "distribution", "window"]
CONSTRUCT_CLASSES = {"wt": "wild type", "mut": "stem mutant"}
MEASUREMENTS = {"frac_1edit": "molecules with exactly 1 edit",
                "frac_2edit": "molecules with 2 or more edits",
                "frac_site": "molecules edited at the recorder site"}
LOOP_MEASUREMENTS = {"frac_1edit": "loop molecules with exactly 1 edit",
                     "frac_2edit": "loop molecules with 2 or 3 edits"}
# How to read the count columns of a panel whose measurement averages two recorder sites: the
# percentage is the geometric mean of the sites' own ratios, so the pooled counts bound it rather
# than reproduce it. The heat-map panels carry each site's own counts beside the average as well.
SITE_BASIS = "counts at the named recorder site"
MEAN_BASIS = "the smaller of the A8 and A2 counts; the plotted percentage is the geometric mean of the two sites"
CONTROL_BASIS = "a control's ratio is to the pooled wild type of its own site, then averaged over the sites"
MUTANT_SERIES = "designed sequence"
CONTROL_SERIES = "wild-type control"
# The pooled matched-wild-type effect columns, in the words a reader outside this analysis needs.
EFFECT_NAMES = {"pooled_edited_umis": "edited_umis", "pooled_umis": "callable_umis",
                "pooled_matched_wt_edited_umis": "wildtype_edited_umis", "pooled_matched_wt_umis": "wildtype_callable_umis",
                "percent_of_matched_wt": "percent_of_wildtype", "log2_pooled_ratio_to_matched_wt": "log2_ratio_to_wildtype",
                "n_ratio_geometries": "n_placements", "n_ratio_windows": "n_windows", "n_ratio_constructs": "n_constructs"}
EFFECT_COLUMNS = ["edited_umis", "callable_umis", "wildtype_edited_umis", "wildtype_callable_umis",
                  "percent_of_wildtype", "log2_ratio_to_wildtype"]
REACTION_COLUMNS = ["enzyme", "tada_concentration", "incubation"]

labels, _, _ = palettes()
# Each panel names the figure it is the source data for, so a renamed figure fails here rather than
# leaving a table matched to nothing.
panels = annotation("source_data_panels").fillna({"selection": ""})
require(not panels.panel.duplicated().any() and panels.family.isin(FAMILIES).all(), "Invalid source-data panel selections")
require(set(itertools.chain.from_iterable(panels.figure.str.split(";"))) <= set(annotation("figure_outputs").filename),
        "A source-data panel names an unregistered figure")
selections = {family: dict(zip(frame.panel, frame.selection)) for family, frame in panels.groupby("family")}
require(set(selections) == set(FAMILIES), "Every panel family must draw at least one panel")
# The averaged measurement names itself the way its axes do, not by its annotation identifier.
MEAN_LABEL = selected_measurement(MEAN_MEASUREMENT).iloc[0].label
os.makedirs(SOURCE_DIRECTORY, exist_ok=True)
written = []


def save(panel, frame):
    """Write one panel's values, failing rather than shipping an empty or unnamed table."""
    require(len(frame) > 0 and not frame.columns.duplicated().any() and frame.notna().all().all(),
            f"Empty, duplicated, or incomplete source data for {panel}")
    path = f"{SOURCE_DIRECTORY}/{panel}.csv"
    frame.to_csv(path, index=False)
    written.append(path)
    print(f"Saved {path}; {len(frame)} rows over {len(frame.columns)} columns")


# The palette names the enzyme-free reaction after the component it omits, which is how the
# component-omission panels label it. On its own it is simply a reaction carrying no enzyme.
ENZYME_NAMES = {**labels, CONTROL_ENZYME: "No enzyme"}


def reaction_columns(frame):
    """Name the reaction of every row the way every figure of this analysis names it."""
    frame["enzyme"] = frame.tada_type.map(ENZYME_NAMES)
    require(frame.enzyme.notna().all(), "A reaction carries no palette name")
    return frame.rename(columns={"tada_conc": "tada_concentration", "condition": "incubation"})


def effect_columns(frame, measurement, basis):
    """Rename one pooled matched-wild-type effect frame and say what its counts are."""
    frame = reaction_columns(frame.copy()).rename(columns=EFFECT_NAMES)
    frame["measurement"] = measurement
    frame["counts_basis"] = basis
    return frame


# Bar panels. Both measurement modes come from one pass over the recorder counts, and a figure that
# stacks rows takes each row's own definition, named in the column that row varies.
definitions, prepared = comparison_data()
bar_panels = {panel: (family, group) for family in ["edit_count", "site"] for panel, group in selections[family].items()}
for panel, (mode, group) in bar_panels.items():
    rows = [row for row in definitions if row["figure_group"] == group and row["measurement_mode"] == mode]
    require(len(rows) > 0, f"No comparison definition draws {panel}")
    frames = []
    for definition in rows:
        summary = prepared[definition["name"]][0].copy()
        column = definition["column"]
        summary["panel_row"] = definition["row_label"] or definition["name"]
        summary["condition_variable"] = column
        summary["condition_level"] = summary[column]
        summary["condition_label"] = summary[column].map(dict(zip(definition["levels"].split(";"), definition["labels"].split(";"))))
        # A comparison that varies the incubation carries it as its level; the rest carry the one
        # incubation their definition selects, so every file names it the same way.
        summary["incubation"] = summary[column] if column == "condition" else definition["conditions"]
        summary = summary.drop(columns="condition", errors="ignore")
        summary["recorder_site"] = definition.get("position_id", "")
        summary["recorder_context"] = str(definition.get("dna_context", "")).replace("T", "U")
        frames.append(summary)
    data = reaction_columns(pd.concat(frames, ignore_index=True))
    data["measurement"] = data.edit_type.map(MEASUREMENTS)
    data["construct_class"] = data.insert_type.map(CONSTRUCT_CLASSES)
    data = data.rename(columns={"n": "n_constructs", "edited_umis": "total_edited_umis", "callable_umis": "total_callable_umis",
                                "mean": "mean_percent_edited", "se": "sem_percent_edited"})
    # Only a figure that stacks several named rows needs to say which row a bar belongs to.
    stacked = ["panel_row"] if len(rows) > 1 else []
    save(panel, data[[*stacked, *REACTION_COLUMNS, "condition_variable", "condition_level", "condition_label",
                      "recorder_site", "recorder_context", "measurement", "construct_class", "n_constructs",
                      "total_edited_umis", "total_callable_umis", "mean_percent_edited", "sem_percent_edited"]])

# Per-adenosine boxB panel, using the same summaries and retained windows as the figure.
for panel, selected_figure in selections["hairpin"].items():
    reactions = hairpin_site_figures()
    data, positions, constructs = hairpin_site_data(reactions.loc[reactions.figure.eq(selected_figure)])
    data = reaction_columns(data)
    data["construct_class"] = data.insert_type.map(CONSTRUCT_CLASSES)
    save(panel, data)

# Direct-versus-split panel: one row per paired variant and recruitment comparison.
for panel, selected_figure in selections["recruitment"].items():
    selection = annotation("boxb_recruitment_comparison")
    require(selection.source.eq(selected_figure).all(), "Unexpected recruitment panel")
    ids = set(selection.x_sample_id) | set(selection.y_sample_id)
    reactions = annotation("sample_info")
    samples = reactions.loc[reactions.sample_id.isin(ids), SAMPLE_KEYS[:4]].copy()
    samples["enzyme_label"] = samples.tada_type.map(labels)
    samples["panel_order"] = range(len(samples))
    effects, _ = sequence_effects(samples, selected_measurement(selection.measurement.iloc[0]))
    effects = effects.loc[effects.n_mutations.gt(0)]
    records = []
    for row in selection.itertuples():
        x = effects.loc[effects.sample_id.eq(row.x_sample_id)].set_index("hairpin")
        y = effects.loc[effects.sample_id.eq(row.y_sample_id)].set_index("hairpin")
        values = pd.DataFrame({"split_sample": row.x_sample_id, "direct_sample": row.y_sample_id,
                               "split_enzyme": x.enzyme_label,
                               "split_log2_ratio_to_wildtype": x.log2_pooled_ratio_to_matched_wt,
                               "direct_log2_ratio_to_wildtype": y.log2_pooled_ratio_to_matched_wt})
        values = values.loc[np.isfinite(values.select_dtypes(include="number")).all(axis=1)]
        values["split_percent_of_wildtype"] = 100 * 2 ** values.split_log2_ratio_to_wildtype
        values["direct_percent_of_wildtype"] = 100 * 2 ** values.direct_log2_ratio_to_wildtype
        values["measurement"] = MEAN_LABEL
        records.append(values.reset_index())
    save(panel, pd.concat(records, ignore_index=True))

# Recorder placements: every drawn point, both enzymes of one concentration in one file.
distance_summary = distance_data()[0]
for panel, concentration in selections["distance"].items():
    data = reaction_columns(distance_summary.loc[distance_summary.tada_conc.eq(concentration)].copy())
    data["construct_class"] = data.insert_type.map(CONSTRUCT_CLASSES)
    data = data.rename(columns={"position_id": "recorder_site", "context": "recorder_context", "target_dist": "boxb_recorder_distance_nt",
                                "n": "n_constructs", "edited_umis": "total_edited_umis", "callable_umis": "total_callable_umis",
                                "mean_percent": "mean_percent_edited", "se_percent": "sem_percent_edited"})
    save(panel, data[["enzyme", "tada_concentration", "recorder_site", "recorder_context", "boxb_recorder_distance_nt",
                      "construct_class", "n_constructs", "total_edited_umis", "total_callable_umis",
                      "mean_percent_edited", "sem_percent_edited"]].sort_values(
                      ["enzyme", "recorder_site", "boxb_recorder_distance_nt", "construct_class"]))

# Flanking-base context tiles: one pooled percentage per site, enzyme, and context.
context_summary = context_data()[0]
for panel, sites in selections["context"].items():
    data = reaction_columns(context_summary.loc[context_summary.position_id.isin(sites.split(";"))].copy())
    data["percent_edited_adenosine"] = 100 * data.pooled_fraction_central_a_edited
    data = data.rename(columns={"position_id": "recorder_site", "fiveprime": "fiveprime_base", "threeprime": "threeprime_base",
                                "n_site_variants": "n_sequence_variants", "n_geometries": "n_placements",
                                "total_edited_umis": "edited_umis", "total_umis": "callable_umis"})
    save(panel, data[[*REACTION_COLUMNS, "recorder_site", "rna_context", "fiveprime_base", "threeprime_base",
                      "n_placements", "n_sequence_variants", "edited_umis", "callable_umis", "percent_edited_adenosine"]])

# boxB single substitutions: the heat map tiles and the logo stacks both read this one table, and
# each row carries the two sites' own counts beside the average the panel draws.
single_mutants = {set_id: single_mutant_data(set_id)[0]
                  for set_id in dict.fromkeys([*selections["heatmap"].values(), *selections["substitution"].values()])}
for panel, set_id in selections["heatmap"].items():
    tiles = single_mutants[set_id]
    data = effect_columns(tiles.loc[tiles.position_id.eq(MEAN_MEASUREMENT)], MEAN_LABEL, MEAN_BASIS)
    for site in CORRELATION_SITES:
        rows = tiles.loc[tiles.position_id.eq(site)].rename(columns=EFFECT_NAMES)
        rows = rows[["sample_id", "boxb_position", "nucleotide", *EFFECT_COLUMNS]].rename(
            columns={name: f"{site.lower()}_{name}" for name in EFFECT_COLUMNS})
        data = data.merge(rows, on=["sample_id", "boxb_position", "nucleotide"], validate="one_to_one")
    data = data.rename(columns={"nucleotide": "mutant_base"})
    save(panel, data[["sample_id", *REACTION_COLUMNS, "measurement", "counts_basis", "boxb_position", "wildtype_base",
                      "mutant_base", "is_wildtype", "n_placements", "n_windows", "n_constructs", *EFFECT_COLUMNS,
                      *[f"{site.lower()}_{name}" for site in CORRELATION_SITES for name in EFFECT_COLUMNS]]].sort_values(
                      ["sample_id", "boxb_position", "mutant_base"]))

# The A8-against-A2 scatters, over the 57 single substitutions or over every designed sequence.
designed = {set_id: designed_variant_effects(set_id, CORRELATION_SITES) for set_id in dict.fromkeys(selections["designed"].values())}
for panel, set_id in {**selections["designed"], **selections["substitution"]}.items():
    whole_library = panel in selections["designed"]
    source = designed[set_id] if whole_library else single_mutants[set_id]
    keys = ["sample_id", "hairpin"] + ([] if whole_library else ["boxb_position", "nucleotide"])
    kept = ["n_mutations"] + ([] if whole_library else ["wildtype_base", "is_wildtype"])
    frames = []
    for site in CORRELATION_SITES:
        rows = source.loc[source.position_id.eq(site)].rename(columns=EFFECT_NAMES)
        rows = rows[keys + kept + EFFECT_COLUMNS].rename(columns={name: f"{site.lower()}_{name}" for name in EFFECT_COLUMNS})
        frames.append(rows if not frames else rows.drop(columns=kept))
    reactions = reaction_columns(source[["sample_id", "tada_type", "tada_conc", "condition"]].drop_duplicates())
    data = frames[0].merge(frames[1], on=keys, validate="one_to_one").merge(reactions, on="sample_id", validate="many_to_one")
    # Wild type is the reference every ratio is taken against, so the panel leaves it out; a sequence
    # with no edited molecule at a site has no logarithm and is counted in the figure rather than drawn.
    finite = np.isfinite(data[[f"{site.lower()}_log2_ratio_to_wildtype" for site in CORRELATION_SITES]]).all(axis=1)
    data["drawn_in_panel"] = data.n_mutations.gt(0) & finite
    data["counts_basis"] = SITE_BASIS
    data = data.rename(columns={"hairpin": "boxb_sequence", "nucleotide": "mutant_base", "n_mutations": "n_substitutions_from_wildtype"})
    data["boxb_sequence"] = data.boxb_sequence.str.replace("T", "U")
    identifiers = ["sample_id", *REACTION_COLUMNS, "counts_basis", "boxb_sequence", "n_substitutions_from_wildtype"]
    identifiers += [] if whole_library else ["boxb_position", "wildtype_base", "mutant_base"]
    data = data[identifiers + ["drawn_in_panel"] + [f"{site.lower()}_{name}" for site in CORRELATION_SITES for name in EFFECT_COLUMNS]]
    save(panel, data.replace([np.inf, -np.inf], "").sort_values(["sample_id", "boxb_sequence"]))

# The distribution of designed-sequence effects by distance from wild type, with the wild-type
# controls that hold position zero of every panel.
for panel, set_id in selections["distribution"].items():
    distances, controls, measurement = mutation_distance_data(set_id)
    mutants = effect_columns(distances, measurement.label, MEAN_BASIS)
    mutants["series"] = MUTANT_SERIES
    mutants = mutants.rename(columns={"hairpin": "boxb_sequence", "n_mutations": "n_substitutions_from_wildtype"})
    mutants["boxb_sequence"] = mutants.boxb_sequence.str.replace("T", "U")
    points = reaction_columns(controls.copy())
    points["series"] = CONTROL_SERIES
    points["measurement"] = measurement.label
    points["counts_basis"] = CONTROL_BASIS
    points["n_substitutions_from_wildtype"] = 0
    points["percent_of_wildtype"] = 100 * points.ratio_to_pooled_wt
    points = points.rename(columns={"log2_ratio_to_pooled_wt": "log2_ratio_to_wildtype", "barcode": "construct_barcode"})
    shared = ["series", *REACTION_COLUMNS, "measurement", "counts_basis", "n_substitutions_from_wildtype",
              "percent_of_wildtype", "log2_ratio_to_wildtype"]
    data = pd.concat([points[[*shared, "construct_barcode"]],
                      mutants[[*shared, "boxb_sequence", *EFFECT_COLUMNS[:4]]]], ignore_index=True)
    # A sequence with no edited molecule has no logarithm; the violin leaves it out and counts it.
    # Wild type is drawn as its own control points, so the designed wild-type sequence has no violin.
    data["drawn_in_panel"] = np.isfinite(data.log2_ratio_to_wildtype) & ~(data.series.eq(MUTANT_SERIES) & data.n_substitutions_from_wildtype.eq(0))
    save(panel, data.fillna("").replace([np.inf, -np.inf], ""))

# The six combinatorial window maps of one reaction, one row per drawn cell.
for panel, set_id in selections["window"].items():
    cells, _, sample, measurement = multimutant_window_data(set_id)
    data = effect_columns(cells.assign(tada_conc=sample.tada_conc, condition=sample.condition), measurement.label, MEAN_BASIS)
    data["window"] = "boxB " + data.window_start.astype(str) + "-" + data.window_end.astype(str)
    data = data.rename(columns={"window_sequence": "window_bases", "hairpin": "boxb_sequence", "row_label": "map_row",
                                "column_label": "map_column", "window_umis": "window_callable_umis"})
    for column in ["window_bases", "boxb_sequence", "wildtype_sequence"]:
        data[column] = data[column].str.replace("T", "U")
    save(panel, data[[*REACTION_COLUMNS, "measurement", "counts_basis", "window", "map_row", "map_column",
                      "window_bases", "wildtype_sequence", "boxb_sequence", "n_substitutions", "is_wildtype",
                      "n_placements", "n_windows", "window_callable_umis", "window_edited_umis",
                      *EFFECT_COLUMNS]].replace([np.inf, -np.inf], "").sort_values(["window", "map_row", "map_column"]))

require(set(written) == {f"{SOURCE_DIRECTORY}/{panel}.csv" for panel in panels.panel}, "A declared panel was not written")
print(f"Wrote {len(written)} panel tables to {SOURCE_DIRECTORY}/")
