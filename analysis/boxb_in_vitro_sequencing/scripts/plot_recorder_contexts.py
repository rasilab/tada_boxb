"""Depth-weighted recorder context preference, with no coverage cutoff or correction."""

import math
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from analysis_helpers import annotation, context_data, require, RNA_BASES
from plot_style import plt, palettes, save_pdf, place_shared_ylabel, centre_shared_xlabel

PANEL_SIZE_IN = .82
PANEL_GAP_IN = .52
TITLE_PAD_PT = 2
MARGINS_IN = {"left": .55, "bottom": .31, "top": .40}
FIGURE_WIDTH_IN = 1.80
FIGURE_SIZE = (FIGURE_WIDTH_IN, 2 * PANEL_SIZE_IN + PANEL_GAP_IN + MARGINS_IN["bottom"] + MARGINS_IN["top"])
# Sites that share one output are drawn as a grid of tiles, one column per site and one row per
# enzyme, so the two flanking-base axes are named once for the whole grid instead of once per site.
GRID_GAP_IN = .07
GRID_MARGINS_IN = {"left": .45, "bottom": .33, "top": .42, "right": .22}
GRID_ROW_LABEL_GAP = .05
TITLE = "% Edited Adenosine"
LOW_COLOR = "#F7F7F7"
SCALE_STEP = 10
WHITE_TEXT_THRESHOLD = .65

summary, site_summary = context_data()
labels, _, colors = palettes()
site_order = annotation("recorder_context_selection")[["position_id", "position_order"]].drop_duplicates()
site_order = site_order.set_index("position_id", verify_integrity=True).position_order
outputs = annotation("figure_outputs")
outputs = outputs.loc[outputs.analysis.eq("plot_recorder_contexts")]
# A site drawn twice, or not at all, would otherwise show up only as a missing or doubled tile.
drawn = [site for output in outputs.itertuples() for site in output.position_id.split(";")]
require(len(drawn) == len(set(drawn)) and set(drawn) == set(summary.position_id), "Draw each selected recorder site exactly once")


def draw_tile(ax, rows, scale_max):
    """Draw one 4x4 flanking-base tile, its printed percentages shaded by the row's enzyme."""
    enzyme = rows.iloc[0].tada_type
    values = rows.pivot(index="fiveprime", columns="threeprime", values="pooled_fraction_central_a_edited").loc[list(RNA_BASES), list(RNA_BASES)].to_numpy() * 100
    cmap = LinearSegmentedColormap.from_list("enzyme", [LOW_COLOR, colors[enzyme]])
    ax.pcolormesh(np.arange(5) - .5, np.arange(5) - .5, values, cmap=cmap, vmin=0, vmax=scale_max, rasterized=False)
    for (y, x), value in np.ndenumerate(values):
        ax.text(x, y, f"{value:.0f}", ha="center", va="center", color="white" if value > scale_max * WHITE_TEXT_THRESHOLD else "black")
    ax.set(xticks=range(4), xticklabels=list(RNA_BASES), yticks=range(4), yticklabels=list(RNA_BASES), xlim=(-.5, 3.5), ylim=(3.5, -.5), aspect="equal")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return enzyme


for output in outputs.itertuples():
    sites = output.position_id.split(";")
    if len(sites) == 1:
        site_data = summary.loc[summary.position_id.eq(output.position_id)]
        scale_max = max(SCALE_STEP, math.ceil(100 * site_data.pooled_fraction_central_a_edited.max() / SCALE_STEP) * SCALE_STEP)
        figure = plt.figure(figsize=FIGURE_SIZE)
        axes = [figure.add_axes([MARGINS_IN["left"] / FIGURE_SIZE[0],
                                (MARGINS_IN["bottom"] + row * (PANEL_SIZE_IN + PANEL_GAP_IN)) / FIGURE_SIZE[1],
                                PANEL_SIZE_IN / FIGURE_SIZE[0], PANEL_SIZE_IN / FIGURE_SIZE[1]]) for row in [1, 0]]
        figure.suptitle(f"{TITLE}, {output.position_id}")
        for ax, (_, rows) in zip(axes, site_data.groupby("panel_order", sort=True)):
            enzyme = draw_tile(ax, rows, scale_max)
            ax.set_title(labels[enzyme], pad=TITLE_PAD_PT)
        axes[-1].set_xlabel("3′ flanking base")
        place_shared_ylabel(figure, axes[0], "5′ flanking base")
        save_pdf(figure, output.source, transparent=True, png=True)
        continue

    # One shared colour scale across the grid, so a tile reads against its neighbours rather than
    # against itself; the printed percentages stay the quantitative record either way.
    require(set(sites).issubset(site_order.index), "Unregistered recorder site in a combined context figure")
    columns = sorted(sites, key=lambda site: site_order[site])
    grid_data = summary.loc[summary.position_id.isin(columns)]
    enzyme_rows = sorted(grid_data.panel_order.unique())
    scale_max = max(SCALE_STEP, math.ceil(100 * grid_data.pooled_fraction_central_a_edited.max() / SCALE_STEP) * SCALE_STEP)
    width = GRID_MARGINS_IN["left"] + len(columns) * PANEL_SIZE_IN + (len(columns) - 1) * GRID_GAP_IN + GRID_MARGINS_IN["right"]
    height = GRID_MARGINS_IN["bottom"] + len(enzyme_rows) * PANEL_SIZE_IN + (len(enzyme_rows) - 1) * GRID_GAP_IN + GRID_MARGINS_IN["top"]
    figure = plt.figure(figsize=(width, height))
    figure.suptitle(TITLE)
    left_column = []
    for row, panel_order in enumerate(enzyme_rows):
        for column, site in enumerate(columns):
            rows = grid_data.loc[grid_data.panel_order.eq(panel_order) & grid_data.position_id.eq(site)]
            require(len(rows) == len(RNA_BASES) ** 2, "Incomplete flanking-base tile")
            ax = figure.add_axes([(GRID_MARGINS_IN["left"] + column * (PANEL_SIZE_IN + GRID_GAP_IN)) / width,
                                  (GRID_MARGINS_IN["bottom"] + (len(enzyme_rows) - 1 - row) * (PANEL_SIZE_IN + GRID_GAP_IN)) / height,
                                  PANEL_SIZE_IN / width, PANEL_SIZE_IN / height])
            enzyme = draw_tile(ax, rows, scale_max)
            # Name each site once above its column and each enzyme once beside its row, and keep the
            # base ticks only on the outer edges, which is what tightens the grid over separate figures.
            if row == 0:
                ax.set_title(site, pad=TITLE_PAD_PT)
            if row < len(enzyme_rows) - 1:
                ax.set_xticklabels([])
            if column == 0:
                left_column.append(ax)
            else:
                ax.set_yticklabels([])
            if column == len(columns) - 1:
                ax.text(1 + GRID_ROW_LABEL_GAP, .5, labels[enzyme], transform=ax.transAxes, rotation=270, ha="left", va="center")
    place_shared_ylabel(figure, left_column[0], "5′ flanking base")
    centre_shared_xlabel(figure, figure.supxlabel("3′ flanking base"))
    require(all(text.get_window_extent(figure.canvas.get_renderer()).x1 <= figure.bbox.x1
                for ax in figure.axes for text in ax.texts if text.get_rotation() == 270), "An enzyme label does not fit the page")
    save_pdf(figure, output.source, transparent=True, png=True)
