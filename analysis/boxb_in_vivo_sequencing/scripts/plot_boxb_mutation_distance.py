"""Mean recorder editing by substitution count, with an in vitro-style violin layout.

Each mutant contributes its pooled editing percentage averaged equally over the selected sites.
Position zero shows distinct WT constructs. No reference line is drawn.
Sequences and WT constructs are not biological replicates. Observed zero editing is retained.
"""

import os

import numpy as np
import pandas as pd

from analysis_helpers import (annotation, comparison_samples, designed_boxb_sequences, excluded_constructs,
                              require, sequence_effects, substitution_measurements)
from plot_style import palettes, plt, save_figure

FIGURE_SIZE = (1.50, 1.24)
VIOLIN_WIDTH = .80
VIOLIN_ALPHA = .55
MEDIAN_WIDTH = .26
MEDIAN_LINEWIDTH = .9
CONTROL_POINT_SIZE = 1.6
CONTROL_JITTER = .20
CONTROL_SEED = 0
Y_TICK_STEP = 5
X_PADDING = .65
TABLE_DIRECTORY = "analysis/boxb_in_vivo_sequencing/data/figure_data"

_, _, colors = palettes()
selections = annotation("boxb_mutation_distance")
require(selections.notna().all().all() and not selections.source.duplicated().any(), "Invalid violin selections")
require(excluded_constructs().empty, "Flagged placements require an explicit with/without comparison")
designed = designed_boxb_sequences()

for selection in selections.itertuples():
    selected = comparison_samples(selection.comparison_set)
    selected = selected.loc[selected.sample_id.eq(selection.sample_id)]
    require(len(selected) == 1, "A violin panel must select one annotated library")
    sites = selection.sites.split(";")
    measurements = substitution_measurements(selection.comparison_set).set_index("position_id")
    require(len(sites) == len(set(sites)) and set(sites) <= set(measurements.index), "Invalid violin sites")
    measurements = measurements.loc[sites].reset_index()
    require(measurements.sites.eq(measurements.position_id).all(), "Average marginal sites, not joint events")
    effects, constructs = sequence_effects(selected, measurements)
    effects = effects.merge(designed[["hairpin", "n_mutations"]], on=["hairpin", "n_mutations"], validate="many_to_one")
    require(len(effects) == len(designed) * len(sites)
            and not effects.duplicated(["hairpin", "position_id"]).any(), "Missing designed sequence measurement")
    require(effects.pooled_umis.gt(0).all() and effects.pooled_fraction_edited.between(0, 1).all(),
            "A designed sequence lacks callable coverage or valid editing")
    average = effects.groupby(["sample_id", "hairpin", "n_mutations"], as_index=False).agg(
        mean_fraction_edited=("pooled_fraction_edited", "mean"), n_sites=("position_id", "size"))
    require(average.n_sites.eq(len(sites)).all(), "Every sequence must include every selected site")
    average["mean_editing_percent"] = 100 * average.mean_fraction_edited

    controls = constructs.loc[constructs.n_mutations.eq(0) & constructs.included_in_ratio].copy()
    require(not controls.duplicated(["barcode", "position_id"]).any()
            and controls.callable_umis.gt(0).all(), "Invalid WT construct measurements")
    controls["mean_editing_percent"] = 100 * controls.edited_umis / controls.callable_umis
    controls = controls.groupby(["sample_id", "hairpin", "barcode", "n_mutations"], as_index=False).agg(
        mean_editing_percent=("mean_editing_percent", "mean"), n_sites=("position_id", "size"))
    require(controls.n_sites.eq(len(sites)).all(), "Every WT construct must include every selected site")

    distances = sorted(average.n_mutations.unique())
    violins = [controls.mean_editing_percent.to_numpy()] + [
        average.loc[average.n_mutations.eq(distance), "mean_editing_percent"].to_numpy()
        for distance in distances if distance > 0]
    require(distances[0] == 0 and all(len(values) > 1 and np.ptp(values) > 0 for values in violins),
            "A substitution count has insufficient variation to form a violin")
    library = selected.iloc[0]
    color = colors[library.enzyme]
    figure, ax = plt.subplots(figsize=FIGURE_SIZE, layout="constrained")
    pooled_wt = average.loc[average.n_mutations.eq(0), "mean_editing_percent"].item()
    parts = ax.violinplot(violins, positions=distances, widths=VIOLIN_WIDTH, showextrema=False, showmedians=False)
    for body in parts["bodies"]:
        body.set(facecolor=color, edgecolor="none", alpha=VIOLIN_ALPHA)
    offsets = np.random.default_rng(CONTROL_SEED).uniform(-CONTROL_JITTER, CONTROL_JITTER, len(controls))
    ax.plot(offsets, controls.mean_editing_percent, linestyle="none", marker="o", markersize=CONTROL_POINT_SIZE,
            markerfacecolor=color, markeredgecolor="none")
    for distance, values in zip(distances, violins):
        median = np.median(values)
        ax.plot([distance - MEDIAN_WIDTH, distance + MEDIAN_WIDTH], [median, median], color=color,
                linewidth=MEDIAN_LINEWIDTH, solid_capstyle="butt")
        print(f"{distance} substitutions: n={len(values)}, median={median:.3f}%, "
              f"IQR={np.percentile(values, 25):.3f}-{np.percentile(values, 75):.3f}%, "
              f"observed zeros={int((values == 0).sum())}")
    upper = np.ceil(np.concatenate(violins).max() / Y_TICK_STEP) * Y_TICK_STEP
    ticks = np.arange(0, upper + Y_TICK_STEP, Y_TICK_STEP)
    ax.set(xlabel="Substitutions from\nwild type",
           ylabel=selection.label.replace(" editing", "\nediting"),
           xticks=distances, xlim=(-X_PADDING, max(distances) + X_PADDING),
           ylim=(0, ticks[-1]), yticks=ticks)
    ax.set_box_aspect(1)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    for text in [ax.xaxis.label, ax.yaxis.label, ax.title, *ax.get_xticklabels(), *ax.get_yticklabels()]:
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            require(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width
                    and box.y1 <= figure.bbox.height, "A violin label extends outside the page")
    os.makedirs(TABLE_DIRECTORY, exist_ok=True)
    pd.concat([average.assign(measurement="pooled_sequence"), controls.assign(measurement="WT_construct")],
              ignore_index=True).to_csv(f"{TABLE_DIRECTORY}/{selection.source}.csv", index=False)
    print(f"Pooled WT reference: {pooled_wt:.3f}%; no placements excluded; with/without estimates identical.")
    save_figure(figure, selection.source)
