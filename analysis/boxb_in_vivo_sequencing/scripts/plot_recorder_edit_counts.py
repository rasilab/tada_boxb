"""Registered WT/MUT edit-count bars with matched-construct tests, one figure per annotated panel.

Bars are distinct-construct means and SEM, not biological replicates. Fold labels divide the
displayed WT mean by the mutant mean. Every denominator requires all eight recorder sites callable.
Paired tests average mutants sharing a barcode before comparing them with their WT, avoiding repeated
WT observations. Tests are uncorrected and describe construct variation.
"""

import os

import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from scipy.stats import ttest_rel

from analysis_helpers import (annotation, comparison_samples, edit_count_fractions, excluded_constructs,
                              matched_stem_pairs, require, substitution_measurements,
                              valid_counts)
from in_vitro_reuse import helpers as in_vitro_helpers
from plot_style import FONT_SIZE, MUTANT_COLORS, palettes, plt, save_figure

EDIT_BINS = {"frac_1edit": "1 edit", "frac_2edit": "2+ edits"}
COUNT_COLUMNS = [f"num_{number}_C" for number in range(1, 9)]
CONSTRUCT_KEYS = ["sample_id", "barcode", "insert"]
GROUP_SPACING = .96
PAIR_OFFSET = .20
BAR_WIDTH = .22
BAR_EDGE_WIDTH = .2
ERROR_LINE_WIDTH = .35
LINE_COLOR = "black"
X_MARGIN = .47
ERROR_CAPSIZE = 1.5
FOLD_DECIMALS = 1
FOLD_GAP_PT = .5
TOP_PADDING_PT = 2
GROUP_LABEL_PAD_PT = 14
LABEL_CLEARANCE_PT = .5
Y_TICK_BINS = 3
SIGNIFICANCE_THRESHOLDS = [.001, .01, .05]
SIGNIFICANCE_LABELS = ["***", "**", "*"]
TABLE_DIRECTORY = "analysis/boxb_in_vivo_sequencing/data/figure_data"
labels, light, dark = palettes()
panels = annotation("recorder_edit_count_panels")
require(not panels.source.duplicated().any(), "Repeated edit-count figure selection")

for panel in panels.itertuples():
    selected = comparison_samples(panel.comparison_set)
    sample_ids = panel.sample_ids.split(";")
    require(len(sample_ids) == len(set(sample_ids)) and set(sample_ids) <= set(selected.sample_id),
            "Repeated or unannotated edit-count library")
    selected = selected.set_index("sample_id").loc[sample_ids].reset_index()
    require(len(selected) > 0 and selected.batch.nunique() == 1, "Missing libraries or mixed batches")
    pairs = matched_stem_pairs(selected, substitution_measurements(panel.comparison_set).iloc[:1])
    bins = edit_count_fractions(selected)
    valid_counts(bins, COUNT_COLUMNS, positive_coverage=False)
    constructs = pairs.drop_duplicates(CONSTRUCT_KEYS)[CONSTRUCT_KEYS + ["insert_type"]]
    counts = constructs.merge(bins, on=CONSTRUCT_KEYS, how="left", validate="one_to_one")
    require(counts[list(EDIT_BINS)].notna().all().all(), "A selected construct has no fully callable molecule")
    require(excluded_constructs().empty, "Flagged placements require an explicit with/without comparison")
    print("No in vivo placements are excluded; with/without-exclusion estimates are identical.")

    long = counts.melt(id_vars=["sample_id", "enzyme", "panel_order", "insert_type"],
                       value_vars=list(EDIT_BINS), var_name="edit_type", value_name="fraction")
    summary = long.groupby(["sample_id", "enzyme", "panel_order", "edit_type", "insert_type"],
                           as_index=False).fraction.agg(n="size", mean="mean", se="sem")
    summary[["mean", "se"]] *= 100
    require(np.isfinite(summary[["mean", "se"]]).all().all(), "Invalid construct mean or SEM")
    matched = pairs.loc[pairs.mutant_bases.notna(), CONSTRUCT_KEYS + ["insert_type"]].merge(
        bins[CONSTRUCT_KEYS + list(EDIT_BINS)], on=CONSTRUCT_KEYS, validate="many_to_one")
    matched = matched.melt(id_vars=CONSTRUCT_KEYS + ["insert_type"], value_vars=list(EDIT_BINS),
                          var_name="edit_type", value_name="fraction")
    require(matched.groupby(["sample_id", "barcode", "edit_type", "insert_type"]).fraction.count().gt(0).all(),
            "Missing matched construct values")
    paired = matched.groupby(["sample_id", "barcode", "edit_type", "insert_type"]).fraction.mean().unstack("insert_type").reset_index()
    require(paired[["wt", "mut"]].notna().all().all(), "A barcode is missing its matched WT or mutants")
    statistics = []
    for (sample_id, edit_type), group in paired.groupby(["sample_id", "edit_type"]):
        test = ttest_rel(group.wt, group.mut)
        statistics.append(dict(sample_id=sample_id, edit_type=edit_type, n_matched_groups=len(group),
                               p_value=test.pvalue, statistic=test.statistic,
                               matched_wt_mean_percent=100 * group.wt.mean(),
                               matched_mut_mean_percent=100 * group.mut.mean()))
    statistics = pd.DataFrame(statistics)
    require(np.isfinite(statistics.p_value).all(), "A paired test is not finite")
    statistics["significance"] = np.select(
        [statistics.p_value.lt(threshold) for threshold in SIGNIFICANCE_THRESHOLDS],
        SIGNIFICANCE_LABELS, default="ns")
    os.makedirs(TABLE_DIRECTORY, exist_ok=True)
    statistics.to_csv(f"{TABLE_DIRECTORY}/{panel.source}_statistics.csv", index=False)
    print(statistics.to_string(index=False))
    figure, axes = plt.subplots(1, len(selected), figsize=(panel.width, panel.height),
                                layout="constrained", sharey=True, squeeze=False)
    annotations = []
    for ax, sample in zip(axes.flat, selected.itertuples()):
        data = summary.loc[summary.sample_id.eq(sample.sample_id)]
        ticks, tick_labels = [], []
        for index, column in enumerate(EDIT_BINS):
            centre = index * GROUP_SPACING
            values = data.loc[data.edit_type.eq(column)].set_index("insert_type")
            shade = "light" if index == 0 else "dark"
            for offset, variant in [(-PAIR_OFFSET, "wt"), (PAIR_OFFSET, "mut")]:
                row = values.loc[variant]
                color = (light if index == 0 else dark)[sample.enzyme] if variant == "wt" else MUTANT_COLORS[shade]
                ax.bar(centre + offset, row["mean"], yerr=row.se, width=BAR_WIDTH,
                       color=color, edgecolor=LINE_COLOR, linewidth=BAR_EDGE_WIDTH,
                       error_kw=dict(elinewidth=ERROR_LINE_WIDTH, capsize=ERROR_CAPSIZE, capthick=ERROR_LINE_WIDTH))
                ticks.append(centre + offset)
                tick_labels.append(variant.upper())
            # A bin whose mutants carry no edits has no ratio to label; the printed summary says so.
            fold = values.loc["wt", "mean"] / values.loc["mut", "mean"] if values.loc["mut", "mean"] > 0 else np.nan
            top = (values["mean"] + values.se).max()
            test = statistics.loc[statistics.sample_id.eq(sample.sample_id) & statistics.edit_type.eq(column)].iloc[0]
            fold_offset = panel.star_offset
            if sample.sample_id not in panel.fold_labels.split(";"):
                labelled = False
            elif test.significance != "ns":
                annotations.append(ax.annotate(test.significance, (centre, top), xytext=(0, panel.star_offset),
                                               textcoords="offset points", ha="center", va="bottom"))
                fold_offset += FONT_SIZE + FOLD_GAP_PT
                labelled = True
            else:
                labelled = True
            if labelled and np.isfinite(fold):
                annotations.append(ax.annotate(f"{fold:.{FOLD_DECIMALS}f}x", (centre, top),
                                               xytext=(0, fold_offset), textcoords="offset points",
                                               ha="center", va="bottom"))
            print(f"{sample.sample_id} {labels[sample.enzyme]}, {EDIT_BINS[column]}: "
                  f"WT/MUT = {fold:.3f} (undefined when MUT has no edits); " + "; ".join(
                      f"{variant.upper()} {values.loc[variant, 'mean']:.3f} ± {values.loc[variant, 'se']:.3f}% "
                      f"(n={int(values.loc[variant, 'n'])})" for variant in ["wt", "mut"]))
        ax.set(xticks=ticks, xticklabels=tick_labels,
               xlim=(-X_MARGIN, GROUP_SPACING + X_MARGIN), ylim=(0, (summary["mean"] + summary.se).max()))
        ax.set_xticks(np.arange(len(EDIT_BINS)) * GROUP_SPACING, list(EDIT_BINS.values()), minor=True)
        ax.tick_params(axis="x", which="major", length=0)
        ax.tick_params(axis="x", which="minor", length=0, pad=GROUP_LABEL_PAD_PT)
        # A figure drawing more than one library names each panel for its enzyme; a figure drawing one
        # names it in the filename alone, as the rest of this analysis does.
        if len(selected) > 1:
            ax.set_title(labels[sample.enzyme])
        if ax is not axes[0, 0]:
            ax.tick_params(axis="y", left=False)
            ax.spines["left"].set_visible(False)
    axes[0, 0].set_ylabel("% Edited RNA")

    # Fit the common y scale to the actual annotation height and verify horizontal label clearance.
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    upper = axes[0, 0].get_ylim()[1]
    for text in annotations:
        height = text.axes.get_window_extent(renderer).height
        reserve = (text.xyann[1] + TOP_PADDING_PT) * figure.dpi / 72 + text.get_window_extent(renderer).height
        require(height > reserve, "Figure is too short for the fold labels")
        upper = max(upper, text.xy[1] / (1 - reserve / height))
    y_ticks = MaxNLocator(nbins=Y_TICK_BINS).tick_values(0, upper)
    axes[0, 0].set_ylim(0, upper)
    axes[0, 0].set_yticks(y_ticks[y_ticks <= upper])
    # A fitted limit can introduce a wider tick label; resolve the margins again before freezing.
    figure.canvas.draw()
    figure.set_layout_engine("none")
    renderer = figure.canvas.get_renderer()
    for ax in axes.flat:
        boxes = [text.get_window_extent(renderer) for text in ax.get_xticklabels()]
        require(all(right.x0 - left.x1 >= LABEL_CLEARANCE_PT * figure.dpi / 72
                    for left, right in zip(boxes, boxes[1:])), "WT/MUT labels overlap")
        for text in [ax.yaxis.label, ax.title, *ax.get_xticklabels(), *ax.get_xticklabels(minor=True),
                     *ax.get_yticklabels(), *ax.texts]:
            if text.get_text() and text.get_visible():
                box = text.get_window_extent(renderer)
                require(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width
                            and box.y1 <= figure.bbox.height,
                            f"Figure label {text.get_text()!r} extends outside the page: {box.bounds}")
    save_figure(figure, panel.source)
