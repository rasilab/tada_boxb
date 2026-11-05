"""BoxB adenosines: WT/MUT editing at each one, read against a reaction carrying no TadA.

The recorder-placement figure's layout read along the hairpin instead of along the recorder, over every
adenosine it holds: the three in the loop and the one that closes the stem. The concentration series is
one grid, a column per editor and a row per concentration, with the no-enzyme reaction the two columns
share centred above them. pAG-TadA and TadA-GFPNb keep a row each, since one carries its own matched
control and the other a single concentration. Every panel takes its own y scale and its own tick
numbers, so a 20-fold concentration range does not flatten the weakest panel.

Every figure opens on a reaction carrying no TadA, because that panel is what the rest has to be read
against. The wild-type/mutant contrast at the two AAG-context adenosines, A10 and A13, is reproduced in
full with no enzyme present, so the contrast alone is not evidence of editing: only its excess over
that first panel is. The four recorder placements are pooled as replicates, so each position carries one
WT and one MUT point; the 5'/3' contrasts that say how far pooling them is fair stay in the statistics
table beside the WT/MUT tests. The WT/MUT stars come from an unpaired two-sided Welch t-test over
exactly the constructs each panel draws, so they test the plotted means and construct counts, and
are reported without correction.
"""

import os
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from analysis_helpers import (hairpin_site_data, hairpin_site_figures, recorder_site_tests, require,
                              CONTROL_ENZYME)
from plot_style import plt, palettes, place_shared_ylabel, LABEL_GAP, MUTANT_COLORS, save_pdf

PANEL_SIZE_IN = (.95, .68)
MARGINS_IN = {"left": .50, "right": .22, "top": .32, "bottom": .34}
COLUMN_GAP_IN = .30
ROW_GAP_IN = .18
# Extra room between a grid's control row and the editor columns under it, for the column headings
# and for the control's own position rows, which it carries because nothing lines up beneath it. The
# control reads as its own measurement rather than the first row of the grid, so the gap is generous.
HEADING_GAP_IN = .46
INSERT_OFFSETS = {"wt": -.18, "mut": .18}
# A panel's wild-type points carry its editor's colour, so a grid reads as one block per editor; the
# shared no-TadA panel belongs to neither and takes black, which is also what the legend names.
CONTROL_COLOR = "black"
Y_LABEL = "% Edited Adenosine"
ANNOTATION_ROWS_PT = {"A context": -10, "A position": -18.5}
ROW_LABEL_PAD_PT = 2
SIDE_LABEL_PAD_PT = 3
LEGEND_HANDLE_PAD_PT = .2
LEGEND_COLUMN_SPACING = 1.0
MARKER_SIZE = 2.5
ERROR_WIDTH = .4
ERROR_CAPSIZE = 1
BRACKET_GAP_PT = 3
BRACKET_LEG_PT = 2
LABEL_HEIGHT_PT = 9
STAR_OFFSET_PT = -2
BRACKET_WIDTH = .4
TICK_BINS = 4
TABLE_DIRECTORY = "analysis/boxb_in_vitro_sequencing/data/figure_data"
TESTS_PER_POSITION = 3

labels, _, dark_colors = palettes()
os.makedirs(TABLE_DIRECTORY, exist_ok=True)

for name, reactions in hairpin_site_figures().groupby("figure", sort=False):
    summary_data, positions, constructs = hairpin_site_data(reactions)
    statistics = recorder_site_tests(constructs, plotted=True)
    require(statistics.groupby("sample_id").size().eq(TESTS_PER_POSITION * len(positions)).all(),
            f"Incomplete hairpin site comparisons in {name}")
    # The stars test the constructs the points draw: same counts and same means, WT and MUT alike.
    drawn = summary_data.pivot(index=["sample_id", "position_id"], columns="insert_type", values=["n", "mean_percent"])
    drawn.columns = [f"{value}_{variant}" for value, variant in drawn.columns]
    tested = statistics.loc[statistics.comparison.eq("WT/MUT")].merge(
        drawn, left_on=["sample_id", "position_id"], right_index=True, validate="one_to_one")
    require(tested.n_first.eq(tested.n_wt).all() and tested.n_second.eq(tested.n_mut).all()
            and np.allclose(tested.first_mean_percent, tested.mean_percent_wt)
            and np.allclose(tested.second_mean_percent, tested.mean_percent_mut),
            f"The WT/MUT tests in {name} do not use the plotted constructs")
    panels = reactions.sort_values(["panel_row", "panel_column"])
    row_widths = panels.groupby("panel_row").size()
    columns, grid = row_widths.max(), len(row_widths) > 1
    control = panels.sample_id.iloc[0]
    wt_colors = {enzyme: CONTROL_COLOR if enzyme == CONTROL_ENZYME else dark_colors[enzyme]
                 for enzyme in panels.enzyme}

    # A row holding fewer panels than the widest is centred over the panel area, which is where the
    # no-TadA reaction two editor columns share sits.
    area_width = columns * PANEL_SIZE_IN[0] + (columns - 1) * COLUMN_GAP_IN
    figure_size = (MARGINS_IN["left"] + area_width + MARGINS_IN["right"],
                   MARGINS_IN["top"] + len(row_widths) * PANEL_SIZE_IN[1]
                   + (len(row_widths) - 1) * ROW_GAP_IN + grid * HEADING_GAP_IN + MARGINS_IN["bottom"])
    figure = plt.figure(figsize=figure_size)
    axes = {}
    for panel in panels.itertuples():
        across = row_widths[panel.panel_row]
        row_width = across * PANEL_SIZE_IN[0] + (across - 1) * COLUMN_GAP_IN
        left = (MARGINS_IN["left"] + (area_width - row_width) / 2
                + panel.panel_column * (PANEL_SIZE_IN[0] + COLUMN_GAP_IN))
        top = (figure_size[1] - MARGINS_IN["top"] - panel.panel_row * (PANEL_SIZE_IN[1] + ROW_GAP_IN)
               - (HEADING_GAP_IN if grid and panel.panel_row else 0))
        axes[panel.sample_id] = figure.add_axes([left / figure_size[0], (top - PANEL_SIZE_IN[1]) / figure_size[1],
                                                PANEL_SIZE_IN[0] / figure_size[0],
                                                PANEL_SIZE_IN[1] / figure_size[1]])

    for panel in panels.itertuples():
        ax = axes[panel.sample_id]
        rows = summary_data.loc[summary_data.sample_id.eq(panel.sample_id)]
        for variant, offset in INSERT_OFFSETS.items():
            points = rows.loc[rows.insert_type.eq(variant)].sort_values("position_index")
            require(len(points) == len(positions), f"A panel of {name} is missing a position")
            color = wt_colors[panel.enzyme] if variant == "wt" else MUTANT_COLORS["dark"]
            ax.errorbar(points.position_index + offset, points.mean_percent, yerr=points.se_percent,
                        fmt="o", markersize=MARKER_SIZE, linewidth=ERROR_WIDTH, capsize=ERROR_CAPSIZE, color=color)
        ax.set_xlim(.55, len(positions) + .45)
        ax.set_xticks([])

    # Each panel keeps its own scale and its own tick numbers, so each takes the headroom its own
    # brackets need and the brackets are then drawn against that final scale.
    figure.canvas.draw()
    height_pt = next(iter(axes.values())).get_window_extent().height * 72 / figure.dpi
    level_pt = BRACKET_GAP_PT + BRACKET_LEG_PT + LABEL_HEIGHT_PT
    tests_drawn, stars = [], {}
    for panel in panels.itertuples():
        rows = summary_data.loc[summary_data.sample_id.eq(panel.sample_id)]
        tests = statistics.loc[statistics.sample_id.eq(panel.sample_id)]
        tests_drawn.append(tests.assign(figure=name, panel_label=panel.panel_label, enzyme=panel.enzyme,
                                        tada_type=panel.tada_type, tada_conc=panel.tada_conc))
        tops = rows.assign(top=rows.mean_percent + rows.se_percent).groupby("position_id").top.max()
        marked = tests.loc[tests.stars.ne("ns") & tests.comparison.eq("WT/MUT")].set_index("position_id").stars
        stars[panel.sample_id] = (tops, marked)
        axes[panel.sample_id].set_ylim(0, max([tops.max()]
                                              + [tops[position] / (1 - level_pt / height_pt) for position in marked.index]))

    for panel in panels.itertuples():
        ax = axes[panel.sample_id]
        ax.yaxis.set_major_locator(MaxNLocator(nbins=TICK_BINS))
        per_pt = ax.get_ylim()[1] / height_pt
        tops, marked = stars[panel.sample_id]
        for position in positions.loc[positions.position_id.isin(marked.index)].itertuples():
            edges = [position.position_index + INSERT_OFFSETS[variant] for variant in INSERT_OFFSETS]
            y = tops[position.position_id] + per_pt * (BRACKET_GAP_PT + BRACKET_LEG_PT)
            ax.plot([edges[0], edges[0], edges[1], edges[1]],
                    [y - per_pt * BRACKET_LEG_PT, y, y, y - per_pt * BRACKET_LEG_PT],
                    color="black", linewidth=BRACKET_WIDTH, solid_capstyle="butt")
            ax.annotate(marked[position.position_id], (sum(edges) / 2, y), xytext=(0, STAR_OFFSET_PT),
                        textcoords="offset points", ha="center", va="bottom")

    # One row of panels heads each with its own concentration and names the editor once at its right.
    # A grid is the transpose: the editor heads its column, over the first row that carries one, and
    # the concentration names its row at the right, as the context heat-map grids name their enzymes.
    for panel in panels.itertuples():
        if not grid or panel.panel_row == 0:
            axes[panel.sample_id].set_title(panel.panel_label)
        elif panel.panel_row == 1:
            axes[panel.sample_id].set_title(labels[panel.enzyme])
    for panel in panels.loc[panels.groupby("panel_row").panel_column.idxmax()].itertuples():
        if grid and panel.panel_row == 0:
            continue
        axes[panel.sample_id].annotate(panel.panel_label if grid else labels[panel.enzyme], (1, .5),
                                       xycoords="axes fraction", xytext=(SIDE_LABEL_PAD_PT, 0),
                                       textcoords="offset points", rotation=270, ha="left", va="center")

    # The position rows sit under the bottom panel of each column, as the stacked edit-count figures
    # carry theirs, so a grid names its adenosines once per column rather than once per panel. A row
    # holding fewer panels than the widest is centred over them, so its x positions sit above no column
    # and the rows at the foot of the grid read for nothing: it carries its own instead.
    bottom = panels.loc[panels.groupby("panel_column").panel_row.idxmax()]
    centred = panels.loc[panels.panel_row.map(row_widths).lt(columns)]
    annotated = pd.concat([bottom, centred]).drop_duplicates("sample_id")
    for panel in annotated.itertuples():
        ax = axes[panel.sample_id]
        for position in positions.itertuples():
            for label, value in [("A context", position.context), ("A position", position.position_id[1:])]:
                ax.annotate(value, (position.position_index, 0), xycoords=ax.get_xaxis_transform(),
                            xytext=(0, ANNOTATION_ROWS_PT[label]), textcoords="offset points", ha="center", va="center")
    # The row names hang off the left spine of a panel rather than off a data coordinate, so the column
    # they take is the one the margin has to hold, whatever x limits a panel carries. Each band of
    # position rows names itself, so a centred control row reads without the foot of the grid.
    row_names = []
    for _, band in annotated.groupby("panel_row"):
        corner = axes[band.loc[band.panel_column.idxmin(), "sample_id"]]
        row_names += [corner.annotate(label, (0, 0), xycoords="axes fraction", xytext=(-ROW_LABEL_PAD_PT, offset),
                                      textcoords="offset points", ha="right", va="center", annotation_clip=False)
                      for label, offset in ANNOTATION_ROWS_PT.items()]

    handles = [Line2D([], [], marker="o", linestyle="none", markersize=MARKER_SIZE, label=variant.upper(),
                      color=wt_colors[panels.enzyme.iloc[0]] if variant == "wt" else MUTANT_COLORS["dark"])
               for variant in INSERT_OFFSETS]
    panel_centre = (MARGINS_IN["left"] + area_width / 2) / figure_size[0]
    figure.legend(handles=handles, loc="upper center", bbox_to_anchor=(panel_centre, 1), ncol=len(handles),
                  frameon=False, handletextpad=LEGEND_HANDLE_PAD_PT, columnspacing=LEGEND_COLUMN_SPACING,
                  borderaxespad=0)
    # Every panel of the left column carries its own tick numbers, so the shared label is placed
    # against whichever of them reaches furthest left rather than against the top one.
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    edge = min(ax.get_position().x0 for ax in axes.values())
    labelled = min([ax for ax in axes.values() if ax.get_position().x0 - edge < 1e-9],
                   key=lambda ax: min(number.get_window_extent(renderer).x0
                                      for number in ax.get_yticklabels() if number.get_text()))
    # The label reads on one line beside a stack of panels. A single row of panels is a shorter strip
    # than the label is long, and a label reaching past the panels runs into the annotation row names
    # beneath them, so there it breaks over two lines instead.
    y_label = place_shared_ylabel(figure, labelled, Y_LABEL)
    spans = [ax.get_position() for ax in axes.values()]
    panel_height = (max(span.y1 for span in spans) - min(span.y0 for span in spans)) * figure_size[1]
    if y_label.get_window_extent(renderer).height / figure.dpi > panel_height:
        place_shared_ylabel(figure, labelled, "\n".join(Y_LABEL.rsplit(" ", 1)))

    # A centred control row shares neither the grid's columns nor its tick numbers, so the shared label
    # reads as the grid's and leaves the control unnamed. The control takes its own copy beside its own
    # numbers, at the same gap, broken over two lines to fit one panel's height, and the shared label
    # then centres on the grid alone rather than on a stack the control is only loosely part of.
    if grid:
        control_axis = axes[control]
        control_label = figure.text(.5, sum(control_axis.get_position().intervaly) / 2,
                                    "\n".join(Y_LABEL.rsplit(" ", 1)), rotation="vertical",
                                    ha="center", va="center", multialignment="center")
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        numbers = [number.get_window_extent(renderer) for number in control_axis.get_yticklabels()
                   if number.get_text()]
        shift = (min(box.x0 for box in numbers) - LABEL_GAP * figure.dpi
                 - control_label.get_window_extent(renderer).x1)
        control_label.set_x(control_label.get_position()[0] + shift / figure.bbox.x1)
        grid_spans = [ax.get_position() for sample, ax in axes.items() if sample != control]
        y_label.set_y((min(span.y0 for span in grid_spans) + max(span.y1 for span in grid_spans)) / 2)

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    require(all(name.get_window_extent(renderer).x0 >= 0 for name in row_names),
            f"The annotation row names do not fit the page in {name}")
    # A locator hands back tick numbers past the panel's own limit, which are never drawn; counting
    # them as ink would condemn a layout over labels the reader never sees.
    drawn = [(text.get_text(), text.get_window_extent(renderer)) for ax in axes.values()
             for text in ax.texts + ax.get_xticklabels() + [ax.title]
             + [number for number, value in zip(ax.get_yticklabels(), ax.get_yticks())
                if ax.get_ylim()[0] <= value <= ax.get_ylim()[1]] if text.get_text()]
    overlaps = [(first_text, second_text) for index, (first_text, left) in enumerate(drawn)
                for second_text, right in drawn[index + 1:] if left.overlaps(right)]
    require(not overlaps, f"Overlapping hairpin-site annotations in {name}: {overlaps[:5]}")

    tests_drawn = pd.concat(tests_drawn, ignore_index=True)
    tests_drawn.to_csv(f"{TABLE_DIRECTORY}/{name}_statistics.csv", index=False)
    support = summary_data.groupby("insert_type").n.agg(["min", "max"])
    print(f"\n{name}: {support.loc['wt', 'min']} wild-type and {support.loc['mut', 'min']} stem-disrupting "
          f"constructs per position, pooled over all four recorder placements.")
    print("Unpaired Welch tests over every plotted construct, uncorrected.")
    print(tests_drawn.loc[tests_drawn.comparison.eq("WT/MUT"),
                          ["tada_type", "panel_label", "position_id", "first_mean_percent", "second_mean_percent",
                           "difference_percentage_points", "n_first", "n_second", "p_value",
                           "p_value_construct_pairs", "p_value_barcode_means", "stars"]].to_string(index=False, float_format="%.4g"))
    floor = summary_data.loc[summary_data.sample_id.eq(control)].set_index(["position_id", "insert_type"]).mean_percent
    excess = summary_data.loc[~summary_data.sample_id.eq(control)].copy()
    excess["over_control"] = [value - floor[(position, variant)] for value, position, variant
                              in zip(excess.mean_percent, excess.position_id, excess.insert_type)]
    print("Excess over the no-TadA panel, percentage points; the WT/MUT contrast is only a measurement "
          "where it exceeds what this control already shows:")
    print(excess.pivot_table(index=["tada_type", "panel_label"], columns=["position_id", "insert_type"],
                             values="over_control").to_string(float_format="%.2f"))
    save_pdf(figure, name, png=True)
