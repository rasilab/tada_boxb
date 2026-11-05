"""Observed read coverage of every designed boxB sequence, pooled across windows and placements.

Count each accepted read once, including PCR duplicates; UMI counts are exported separately.
Wild type is one sequence, pooled over its barcode representations just like every mutant.
"""

import os

import numpy as np
import pandas as pd
from matplotlib.ticker import NullLocator

from analysis_helpers import (COMPLEMENT, annotation, barcode_metadata, boxb_wildtype,
                              comparison_samples, counts_for, designed_boxb_sequences,
                              designed_construct_grid, excluded_constructs, require)
from plot_style import LABEL_GAP, plt, save_figure

FIGURE_SIZE = (1.20, 1.155)
PANEL_ASPECT = 1
POINT_SIZE = 1.5
POINT_COLOR = "black"
Y_LABEL_PAD = 1
WT_LABEL_OFFSET = (8 / 3, 0)
X_TICKS = [1, 380, 760]
X_PADDING = .025
Y_TICK_EXPONENTS = np.array([2, 4, 6])
TABLE_DIRECTORY = "analysis/boxb_in_vivo_sequencing/data/figure_data"
COUNT_COLUMNS = ["read_counts", "umi_counts"]
SUPERSCRIPTS = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

selections = annotation("boxb_read_coverage")
require(selections.notna().all().all() and not selections.source.duplicated().any(), "Invalid coverage selections")
require(excluded_constructs().empty, "Flagged placements require an explicit with/without comparison")
_, _, wildtype = boxb_wildtype()

for selection in selections.itertuples():
    selected = comparison_samples(selection.comparison_set)
    selected = selected.loc[selected.sample_id.eq(selection.sample_id)]
    require(len(selected) == 1, "A coverage panel must select one annotated library")
    designed = designed_boxb_sequences(selection.g_depleted)
    barcodes = barcode_metadata(selection.g_depleted)
    grid = designed_construct_grid(selection.g_depleted).merge(
        barcodes[["barcode", "window_start", "window_end"]], on="barcode", validate="many_to_one")
    grid["hairpin"] = [wildtype[:start - 1] + insert.translate(COMPLEMENT)[::-1] + wildtype[end:]
                       for start, end, insert in zip(grid.window_start, grid.window_end, grid["insert"])]
    require(set(grid.hairpin) == set(designed.hairpin), "Construct grid and designed sequences disagree")

    counts = counts_for(selection.sample_id, [])
    require(set(counts.barcode) <= set(barcode_metadata(None).barcode), "Counts contain an unannotated barcode")
    counts = counts.loc[counts.barcode.isin(barcodes.barcode)]
    matrix = counts[COUNT_COLUMNS].to_numpy()
    require(np.isfinite(matrix).all() and (matrix == np.floor(matrix)).all()
            and (matrix >= 0).all() and counts.read_counts.ge(counts.umi_counts).all(), "Invalid read or UMI counts")
    keys = ["barcode", "oligo_name", "insert"]
    require(len(counts.merge(grid[keys], on=keys, validate="one_to_one")) == len(counts),
            "An observed construct does not match its design")
    constructs = grid.merge(counts[keys + COUNT_COLUMNS], on=keys, how="left", validate="one_to_one")
    constructs[COUNT_COLUMNS] = constructs[COUNT_COLUMNS].fillna(0).astype("int64")
    constructs["observed"] = constructs.read_counts.gt(0)
    coverage = constructs.groupby("hairpin", as_index=False).agg(
        read_counts=("read_counts", "sum"), umi_counts=("umi_counts", "sum"),
        n_designed_constructs=("barcode", "size"), n_observed_constructs=("observed", "sum"),
        n_placements=("layout", "nunique"), n_windows=("variable_subpos", "nunique"))
    coverage = designed.merge(coverage, on="hairpin", how="left", validate="one_to_one")
    require(len(coverage) == len(designed) and coverage.notna().all().all(), "Missing designed sequence coverage")
    require(coverage[COUNT_COLUMNS].sum().eq(counts[COUNT_COLUMNS].sum()).all(), "Pooling changed count totals")
    require(coverage.read_counts.gt(0).all(), "A designed sequence has zero reads; revise the log-axis display explicitly")
    coverage = coverage.sort_values(["read_counts", "hairpin"], ascending=[False, True]).reset_index(drop=True)
    coverage.insert(0, "rank", np.arange(1, len(coverage) + 1))
    coverage.insert(1, "sample_id", selection.sample_id)
    coverage.insert(2, "g_depleted", selection.g_depleted)

    figure, ax = plt.subplots(figsize=FIGURE_SIZE, layout="constrained")
    ax.scatter(coverage["rank"], coverage.read_counts, s=POINT_SIZE, color=POINT_COLOR, edgecolors="none")
    wt = coverage.loc[coverage.n_mutations.eq(0)]
    require(len(wt) == 1, "Wild type must appear exactly once")
    ax.annotate("WT", (wt["rank"].item(), wt.read_counts.item()), xytext=WT_LABEL_OFFSET,
                textcoords="offset points", ha="left", va="center")
    y_ticks = 10.0 ** Y_TICK_EXPONENTS
    require(coverage.read_counts.between(y_ticks[0], y_ticks[-1]).all(), "A variant falls outside the coverage axis")
    ax.set(xlabel="boxB variant rank", ylabel="Reads per variant", yscale="log", xticks=X_TICKS,
           xlim=(1 - X_PADDING * len(coverage), len(coverage) * (1 + X_PADDING)),
           ylim=(y_ticks[0], y_ticks[-1]))
    ax.set_yticks(y_ticks, ["10" + str(exponent).translate(SUPERSCRIPTS) for exponent in Y_TICK_EXPONENTS])
    ax.yaxis.set_minor_locator(NullLocator())
    ax.grid(False, which="both")
    ax.xaxis.labelpad = LABEL_GAP * 72
    ax.yaxis.labelpad = Y_LABEL_PAD
    ax.set_box_aspect(PANEL_ASPECT)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    for text in [ax.xaxis.label, ax.yaxis.label, *ax.texts, *ax.get_xticklabels(), *ax.get_yticklabels()]:
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            require(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width
                    and box.y1 <= figure.bbox.height, "A coverage label extends outside the page")

    os.makedirs(TABLE_DIRECTORY, exist_ok=True)
    coverage.to_csv(f"{TABLE_DIRECTORY}/{selection.source}.csv", index=False)
    constructs.insert(0, "sample_id", selection.sample_id)
    constructs.to_csv(f"{TABLE_DIRECTORY}/{selection.source}_constructs.csv", index=False)
    print(f"{selection.sample_id}: {len(coverage)} sequences, including one WT; "
          f"{len(barcodes)} constructs, {coverage.read_counts.sum():,} accepted reads.")
    print(f"Reads per sequence: min={coverage.read_counts.min():,}, median={coverage.read_counts.median():,.0f}, "
          f"max={coverage.read_counts.max():,}; WT={wt.read_counts.item():,}; "
          f"UMIs={coverage.umi_counts.sum():,}; {int(coverage.read_counts.ge(100).sum())} sequences with at least 100 reads.")
    print("No placements excluded; with/without estimates identical. Counts include PCR duplicates and are not editing fractions.")
    save_figure(figure, selection.source)
