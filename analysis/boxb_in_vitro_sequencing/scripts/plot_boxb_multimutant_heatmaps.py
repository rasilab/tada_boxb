"""Combinatorial boxB windows: the shared pooled effect of every designed window sequence."""

import itertools
import math
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle
from matplotlib.legend_handler import HandlerBase
from analysis_helpers import annotation, multimutant_window_data, reaction_set, require, scale_groups, RNA_BASES
from plot_style import plt, save_pdf, FONT_SIZE

CELL_SIZE = .105
WINDOWS_PER_ROW = 6
PANEL_GAP = .20
LEFT_MARGIN = .30
RIGHT_MARGIN = .20
TOP_MARGIN = .06
PANEL_TITLE_HEIGHT = .16
ROW_GAP = .12
AXIS_LABEL_HEIGHT = .15
TICK_PAD = 1
ROTATED_TICK_HEIGHT = {1: .12, 2: .20}
COLORBAR_LENGTH = 1.2
COLORBAR_WIDTH = .07
COLORBAR_BLOCK = .23
COLORBAR_GAP = .04
COLORBAR_LABEL_GAP = .06
COLORBAR_TICK_LENGTH = 2
COLORBAR_TICK_PAD = 1
LEGEND_COLUMNS = 3
LEGEND_BOTTOM = .04
LEGEND_GAP = .16
LEGEND_LABEL_GAP = .4
LEGEND_COLUMN_GAP = .8
MINIMUM_WIDTH = 3.5
COLOR_STEP = .5
DIVERGING_COLORS = ["#2166AC", "#F7F7F7", "#B2182B"]
ZERO_EDIT_COLOR = "#FFFFFF"
MARKER_COLOR = "#000000"
SEPARATOR_COLOR = "#999999"
KEY_BOUNDARY_COLOR = "#777777"
SEPARATOR_WIDTH = .3
GROUP_SIZE = 4
WT_EDGE_WIDTH = .6
CORNER_FRACTION = .34
ZERO_EDIT_MARKER_SIZE = 2.2
ZERO_EDIT_MARKER_WIDTH = .4
BASE_SIZE = 6
COMBINED_SOURCE = "boxb_multimutant_heatmaps"
GRID_FIELDS = ["log2_pooled_ratio_to_matched_wt", "percent_of_matched_wt", "n_substitutions", "n_ratio_geometries", "n_ratio_windows",
               "pooled_umis", "pooled_edited_umis", "window_umis", "window_edited_umis", "window_placements"]
LEGEND_KEY = CELL_SIZE * 72 / FONT_SIZE

colormap = LinearSegmentedColormap.from_list("effect", DIVERGING_COLORS)


class CellKey(HandlerBase):
    """Draw a legend key as one map cell: its boundary, and the corner mark where the key carries one.

    The map outlines only wild-type cells, so the corner key takes a grey boundary instead of a black
    one. Without it the key is a mark floating in white space, with nothing to show the cell it sits in.
    """

    def __init__(self, edgecolor, corner):
        super().__init__()
        self.edgecolor, self.corner = edgecolor, corner

    def create_artists(self, legend, handle, xdescent, ydescent, width, height, fontsize, transform):
        size = min(width, height)
        left, top = -xdescent, -ydescent + (height + size) / 2
        artists = [Rectangle((left, top - size), size, size, fill=False, edgecolor=self.edgecolor, linewidth=WT_EDGE_WIDTH)]
        if self.corner:
            artists.append(Polygon([(left, top), (left + CORNER_FRACTION * size, top), (left, top - CORNER_FRACTION * size)], closed=True, facecolor=MARKER_COLOR, edgecolor="none"))
        for artist in artists:
            artist.set_transform(transform)
        return artists


def axis_label(members):
    """Name the axis for the boxB positions it carries, in the wording the substitution panels use."""
    return ("Positions " if len(members) > 1 else "Position ") + ", ".join(str(member) for member in members)


def window_label(choice):
    return f"boxB {choice.window_start}–{choice.window_end}"


def draw_panel(ax, panel, label_rows, scale):
    """One window map: square cells, markers, group separators, and labels spelled out in full."""
    values = panel["log2_pooled_ratio_to_matched_wt"]
    # A cell with no edited molecule has no logarithm, so it keeps the background rather than a value.
    ax.set_facecolor(ZERO_EDIT_COLOR)
    ax.pcolormesh(np.arange(len(panel["column_labels"]) + 1) - .5, np.arange(len(panel["row_labels"]) + 1) - .5, np.ma.masked_invalid(values), cmap=colormap, norm=scale["normalizer"], rasterized=False)
    for (y, x), edited in np.ndenumerate(panel["pooled_edited_umis"]):
        if panel["n_substitutions"][y, x] == 0:
            ax.add_patch(Rectangle((x - .5, y - .5), 1, 1, fill=False, edgecolor=MARKER_COLOR, linewidth=WT_EDGE_WIDTH))
        elif panel["n_substitutions"][y, x] == 1:
            ax.add_patch(Polygon([(x - .5, y - .5), (x - .5 + CORNER_FRACTION, y - .5), (x - .5, y - .5 + CORNER_FRACTION)], closed=True, facecolor=MARKER_COLOR, edgecolor="none"))
        if edited == 0:
            ax.plot(x, y, marker="x", color=MARKER_COLOR, markersize=ZERO_EDIT_MARKER_SIZE, markeredgewidth=ZERO_EDIT_MARKER_WIDTH)
    for boundary in range(GROUP_SIZE, len(panel["row_labels"]), GROUP_SIZE):
        ax.axhline(boundary - .5, color=SEPARATOR_COLOR, linewidth=SEPARATOR_WIDTH)
    for boundary in range(GROUP_SIZE, len(panel["column_labels"]), GROUP_SIZE):
        ax.axvline(boundary - .5, color=SEPARATOR_COLOR, linewidth=SEPARATOR_WIDTH)
    # Every map carries the same 16 rows, so only the first map of a row spells them out; the rest
    # keep their own axis label, which names the boxB positions those rows carry.
    ax.set(title=window_label(panel["choice"]),
           xticks=range(len(panel["column_labels"])), xticklabels=panel["column_labels"], xlabel=axis_label(panel["columns"]), xlim=(-.5, len(panel["column_labels"]) - .5),
           yticks=range(len(panel["row_labels"])), yticklabels=panel["row_labels"] if label_rows else [], ylim=(len(panel["row_labels"]) - .5, -.5))
    ax.tick_params(length=0, pad=TICK_PAD, labelsize=BASE_SIZE)
    # Both axes read along the sequence, so turn the column labels onto their side rather than
    # dropping letters: a column label is as long as the row label beside it.
    plt.setp(ax.get_xticklabels(), rotation=90)
    # Each map names its rows immediately to its right, where the next map's cells are not in the way.
    ax.yaxis.set_label_position("right")
    ax.set_ylabel(axis_label(panel["rows"]), rotation=270, va="bottom", labelpad=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def render(layout, source, scale):
    """Draw the chosen windows as rows of maps at equal square cell size, with the colorbar lying flat beneath them."""
    chosen = [panel for row in layout for panel in row]
    # Each key is one cell at cell size, carrying exactly the mark the map draws on it.
    handles = [Rectangle((0, 0), 1, 1), Rectangle((0, 0), 1, 1)]
    handlers = {handles[0]: CellKey(MARKER_COLOR, corner=False), handles[1]: CellKey(KEY_BOUNDARY_COLOR, corner=True)}
    labels = ["Wild type", "Single substitution"]
    # Only name the zero-edit marker when a drawn cell actually carries it.
    if any((panel["pooled_edited_umis"] == 0).any() for panel in chosen):
        handles.append(Line2D([], [], linestyle="none", marker="x", color=MARKER_COLOR, markersize=ZERO_EDIT_MARKER_SIZE, markeredgewidth=ZERO_EDIT_MARKER_WIDTH))
        labels.append("No edited molecule")
    row_widths = [LEFT_MARGIN + sum(CELL_SIZE * len(panel["column_labels"]) for panel in row) + PANEL_GAP * (len(row) - 1) for row in layout]
    # The maps take the full page: the bar lies flat beneath them rather than standing in a band at the
    # right, so the width is the maps' own and the bar costs height instead.
    figure_width = max(MINIMUM_WIDTH, max(row_widths) + RIGHT_MARGIN)
    panel_height = CELL_SIZE * len(chosen[0]["row_labels"])
    row_heights = [PANEL_TITLE_HEIGHT + panel_height + ROTATED_TICK_HEIGHT[max(len(panel["column_labels"][0]) for panel in row)] + AXIS_LABEL_HEIGHT for row in layout]
    maps_height = sum(row_heights) + ROW_GAP * (len(layout) - 1)
    # The bar and the keys share one band beneath the maps, side by side, so the page carries one block
    # of white space under the maps rather than a bar band and a legend band stacked.
    figure_height = TOP_MARGIN + maps_height + COLORBAR_GAP + COLORBAR_BLOCK + LEGEND_BOTTOM
    figure = plt.figure(figsize=(figure_width, figure_height))
    maps, top = [], TOP_MARGIN
    for index, (row, row_height) in enumerate(zip(layout, row_heights)):
        top += ROW_GAP if index else 0
        left = LEFT_MARGIN
        for position, panel in enumerate(row):
            width = CELL_SIZE * len(panel["column_labels"])
            ax = figure.add_axes([left / figure_width, 1 - (top + PANEL_TITLE_HEIGHT + panel_height) / figure_height, width / figure_width, panel_height / figure_height])
            draw_panel(ax, panel, position == 0, scale)
            maps.append(ax)
            left += width + PANEL_GAP
        top += row_height
    # Put the title to the left of the bar and align both marker keys with its centre, so only
    # the tick numbers need space below the bar.
    colorbar_top = TOP_MARGIN + maps_height + COLORBAR_GAP
    colorbar_ax = figure.add_axes([LEFT_MARGIN / figure_width, 1 - (colorbar_top + COLORBAR_WIDTH) / figure_height, COLORBAR_LENGTH / figure_width, COLORBAR_WIDTH / figure_height])
    colorbar = figure.colorbar(ScalarMappable(norm=scale["normalizer"], cmap=colormap), cax=colorbar_ax, ticks=scale["ticks"], orientation="horizontal")
    colorbar.ax.set_xticklabels([f"{100 * 2. ** power:g}" for power in scale["ticks"]])
    colorbar_ax.set_xlim(*scale["range"])
    colorbar_ax.tick_params(length=COLORBAR_TICK_LENGTH, pad=COLORBAR_TICK_PAD)
    colorbar_ax.text(-COLORBAR_LABEL_GAP / COLORBAR_LENGTH, .5, f"{scale['label']}:",
                    transform=colorbar_ax.transAxes, ha="right", va="center")
    legend = figure.legend(handles, labels, handler_map=handlers, loc="center left", bbox_to_anchor=(0, 1 - (colorbar_top + COLORBAR_WIDTH / 2) / figure_height),
                           ncols=LEGEND_COLUMNS, frameon=False, handlelength=LEGEND_KEY, handleheight=LEGEND_KEY, borderpad=0,
                           handletextpad=LEGEND_LABEL_GAP, columnspacing=LEGEND_COLUMN_GAP, borderaxespad=0)
    # The keys follow the bar at a fixed gap, then bar and keys move together until the band they fill
    # is centred: the two cannot be placed by arithmetic, since the key text sets its own width.
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    anchor = legend.get_bbox_to_anchor()
    legend.set_bbox_to_anchor(((colorbar_ax.get_tightbbox(renderer).x1 + LEGEND_GAP * figure.dpi) / figure.bbox.x1, anchor.y0 / figure.bbox.y1))
    figure.canvas.draw()
    bar, box = colorbar_ax.get_tightbbox(renderer), legend.get_window_extent(renderer)
    shift = (figure.bbox.x1 - (box.x1 - bar.x0)) / 2 - bar.x0
    colorbar_ax.set_position([colorbar_ax.get_position().x0 + shift / figure.bbox.x1, *colorbar_ax.get_position().bounds[1:]])
    anchor = legend.get_bbox_to_anchor()
    legend.set_bbox_to_anchor(((anchor.x0 + shift) / figure.bbox.x1, anchor.y0 / figure.bbox.y1))
    figure.canvas.draw()
    bar, box = colorbar_ax.get_tightbbox(renderer), legend.get_window_extent(renderer)
    require(box.x1 <= figure.bbox.x1 and box.y0 >= 0, f"The legend does not fit the page in {source}")
    require(bar.x0 >= 0 and bar.x1 <= figure.bbox.x1 and bar.y0 >= 0, f"The colorbar does not fit the page in {source}")
    require(max(bar.y1, box.y1) < min(ax.get_tightbbox(renderer).y0 for ax in maps), f"Legend band and maps overlap in {source}")
    require(box.x0 > bar.x1, f"Colorbar and legend overlap in {source}")
    print(f"{source}: {figure_width:.2f} × {figure_height:.2f} inches, {len(layout)} row(s) of maps")
    save_pdf(figure, source, transparent=True, png=True)


# Build every concentration before drawing any, so all of them share one colour mapping and one
# colorbar. A colour then means the same percentage of wild type in every map, and the
# concentrations can be read side by side. The shared mapping is the widest any concentration
# needs, which on the current data is the 125 nM mapping, and the bar spans every concentration's
# cells, so no map carries a cell past the end of its bar.
# Only the sets that carry the selected enzyme can draw these maps; an editor set that does not
# contain it is simply not part of this family.
selection_choice = annotation("boxb_multimutant_heatmap_selection").iloc[0]
enzyme = selection_choice.tada_type
# These maps share one colour mapping across every set they draw, so they are scoped to the scale
# group they belong to: a set outside it cannot silently rescale the maps already in the paper.
set_ids = [set_id for set_id in scale_groups()[selection_choice.scale_group] if enzyme in reaction_set(set_id).tada_type.tolist()]
require(len(set_ids) > 0, f"No reaction set carries {enzyme} for the combinatorial maps")
datasets = {set_id: multimutant_window_data(set_id) for set_id in set_ids}
finite = {concentration: cells.log2_pooled_ratio_to_matched_wt[np.isfinite(cells.log2_pooled_ratio_to_matched_wt)]
          for concentration, (cells, _, _, _) in datasets.items()}
# Colour runs in log₂ of the ratio to matched wild type, symmetrically about it, so equal fold
# changes up and down take equal colour and the gradient is spaced evenly over the fold range.
limit = max([COLOR_STEP] + [math.ceil(values.abs().max() / COLOR_STEP) * COLOR_STEP for values in finite.values()])
# The bar shows only the span the cells occupy; no cell reaches the far end of a symmetric scale.
# Its ticks sit on powers of two, labelled as the percentage of wild-type editing a reader can act on.
colorbar_range = (min(values.min() for values in finite.values()), max(values.max() for values in finite.values()))
shared_label = f"{next(iter(datasets.values()))[3].axis_label}, % of wild type"
require(all(measurement.axis_label == next(iter(datasets.values()))[3].axis_label for _, _, _, measurement in datasets.values()), "Every concentration must measure the same recorder site")
scale = dict(normalizer=Normalize(vmin=-limit, vmax=limit), range=colorbar_range,
             ticks=list(range(math.ceil(colorbar_range[0]), math.floor(colorbar_range[1]) + 1)),
             label=shared_label)
print(f"Shared log₂ colour mapping ±{limit:.2f}; bar shown over {100 * 2 ** colorbar_range[0]:.1f}–{100 * 2 ** colorbar_range[1]:.1f}% of wild type, "
      f"ticks at {', '.join(f'{100 * 2. ** power:g}' for power in scale['ticks'])}%")
for set_id, (per_set, _, _, _) in datasets.items():
    values = per_set.log2_pooled_ratio_to_matched_wt
    values = values[np.isfinite(values)]
    print(f"  {set_id} cells occupy {100 * 2 ** values.min():.1f}–{100 * 2 ** values.max():.1f}% of wild type")

# Every map set is drawn once per reaction set, which only its filename records.
for set_id in set_ids:
    cells, selection, sample, measurement = datasets[set_id]
    print(f"{sample.enzyme_label}, {sample.tada_conc}, {sample.condition}; {measurement.label} editing")

    panels = []
    for choice in selection.itertuples():
        window = cells.loc[cells.variable_subpos.eq(choice.variable_subpos)]
        rows = [int(position) for position in str(choice.row_positions).split(";")]
        columns = [int(position) for position in str(choice.column_positions).split(";")]
        row_labels = ["".join(bases) for bases in itertools.product(RNA_BASES, repeat=len(rows))]
        column_labels = ["".join(bases) for bases in itertools.product(RNA_BASES, repeat=len(columns))]
        grids = {field: window.pivot(index="row_label", columns="column_label", values=field).reindex(index=row_labels, columns=column_labels).to_numpy() for field in GRID_FIELDS}
        require(all(grid.shape == (len(row_labels), len(column_labels)) for grid in grids.values()), "Incomplete window matrix")
        require(all(np.isfinite(grid).all() for field, grid in grids.items() if field != "log2_pooled_ratio_to_matched_wt"), "Incomplete window matrix")
        require(np.array_equal(np.isfinite(grids["log2_pooled_ratio_to_matched_wt"]), grids["pooled_edited_umis"] > 0), "Only zero-edit cells may be non-finite")
        panels.append(dict(choice=choice, rows=rows, columns=columns, row_labels=row_labels, column_labels=column_labels, **grids))
        values = grids["percent_of_matched_wt"]
        # Support is reported here rather than on the figure, which carries only its legend.
        print(f"{window_label(choice)}: {values.shape[0]}×{values.shape[1]} cells, "
              f"window placements {int(grids['window_placements'].min())}–{int(grids['window_placements'].max())}, "
              f"pooled placements {int(grids['n_ratio_geometries'].min())}–{int(grids['n_ratio_geometries'].max())} over "
              f"{int(grids['n_ratio_windows'].min())}–{int(grids['n_ratio_windows'].max())} windows, "
              f"window UMIs {int(grids['window_umis'].min())}–{int(grids['window_umis'].max())} with median {int(np.median(grids['window_umis']))}, "
              f"pooled UMIs {int(grids['pooled_umis'].min())}–{int(grids['pooled_umis'].max())}, "
              f"WT fraction={window.iloc[0].pooled_matched_wt_fraction_edited:.10f}, singles={int((grids['n_substitutions'] == 1).sum())}, "
              f"zero-edit cells={int((grids['pooled_edited_umis'] == 0).sum())}, percent {values.min():.10f} to {values.max():.10f}")
    print(f"{len(cells)} cells across {len(panels)} windows, {cells.hairpin.nunique()} distinct sequences, "
          f"{int(cells.pooled_edited_umis.eq(0).sum())} zero-edit cells")
    render([panels[start:start + WINDOWS_PER_ROW] for start in range(0, len(panels), WINDOWS_PER_ROW)], f"{COMBINED_SOURCE}_{set_id}", scale)
