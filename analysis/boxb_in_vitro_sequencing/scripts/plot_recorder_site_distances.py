"""Recorder placements: construct means/SEM with the WT/MUT contrast at each site drawn as stars."""

import os
from matplotlib.lines import Line2D
from analysis_helpers import annotation, distance_data, recorder_site_tests, EXPECTED_DISTANCES, require
from plot_style import plt, palettes, distance_palettes, DISTANCE_SHADES, save_pdf

FIGURE_SIZE = (4.35, 1.75)
MARGINS_IN = {"left": .62, "right": .04, "top": .23, "bottom": .44}
INSERT_OFFSETS = {"wt": -.24, "mut": .24}
DISTANCE_OFFSETS = [-.075, -.025, .025, .075]
LEGEND_TITLE = "boxB-recorder distance"
LEGEND_ANCHOR = (.99, 1.08)
LEGEND_CLEARANCE_PT = 2
Y_LABEL = "% Edited Adenosine"
ANNOTATION_ROWS_PT = {"boxB": -9, "A context": -17.5, "A position": -26}
PAIR_LABEL_PAD_PT = 6
ROW_LABEL_X = .48
MARKER_SIZE = 2.5
ERROR_WIDTH = .4
ERROR_CAPSIZE = 1
Y_PADDING = 1.10
TABLE_DIRECTORY = "analysis/boxb_in_vitro_sequencing/data/figure_data"
EXPECTED_TESTS = 24
BRACKET_GAP_PT = 3
BRACKET_LEG_PT = 2
LABEL_HEIGHT_PT = 9
STAR_OFFSET_PT = -2
BRACKET_WIDTH = .4
ITERATIONS = 20
HEADROOM_STEP = 1.05

summary_data, metadata, positions, constructs = distance_data()
labels, _, _ = palettes()
distance_colors = distance_palettes()
outputs = annotation("figure_outputs")
statistics = recorder_site_tests(constructs)
require(statistics.groupby("sample_id").size().eq(EXPECTED_TESTS).all(), "Incomplete site comparisons")
os.makedirs(TABLE_DIRECTORY, exist_ok=True)

for output in outputs.loc[outputs.analysis.eq("plot_recorder_site_distances")].itertuples():
    enzyme = output.enzyme
    rows = summary_data.loc[summary_data.tada_conc.eq(output.concentration) & summary_data.tada_type.eq(enzyme)]
    ids = constructs.loc[constructs.tada_conc.eq(output.concentration) & constructs.tada_type.eq(enzyme), "sample_id"].unique()
    require(len(ids) == 1, "Distance figure needs exactly one reaction")
    tests = statistics.loc[statistics.sample_id.eq(ids[0])]
    tests.to_csv(f"{TABLE_DIRECTORY}/{output.source}_statistics.csv", index=False)
    figure, ax = plt.subplots(figsize=FIGURE_SIZE)
    figure.subplots_adjust(left=MARGINS_IN["left"] / FIGURE_SIZE[0], right=1 - MARGINS_IN["right"] / FIGURE_SIZE[0],
                           top=1 - MARGINS_IN["top"] / FIGURE_SIZE[1], bottom=MARGINS_IN["bottom"] / FIGURE_SIZE[1])
    for variant, offset in INSERT_OFFSETS.items():
        for i, distance in enumerate(EXPECTED_DISTANCES):
            points = rows.loc[rows.insert_type.eq(variant) & rows.target_dist.eq(distance)].sort_values("position_index")
            color = DISTANCE_SHADES[i] if variant == "mut" else distance_colors[enzyme][i]
            ax.errorbar(points.position_index + offset + DISTANCE_OFFSETS[i], points.mean_percent,
                        yerr=points.se_percent, fmt="o", markersize=MARKER_SIZE, linewidth=ERROR_WIDTH,
                        capsize=ERROR_CAPSIZE, color=color)
    ax.set(title=labels[enzyme], ylabel=Y_LABEL, xlim=(.55, len(positions) + .45))
    ax.tick_params(axis="x", length=0, pad=PAIR_LABEL_PAD_PT)
    handles = [Line2D([], [], marker="o", linestyle="none", color=color, markersize=MARKER_SIZE,
                      label=f"{distance} nt") for color, distance in zip(DISTANCE_SHADES, EXPECTED_DISTANCES)]
    legend = ax.legend(handles=handles, title=LEGEND_TITLE, loc="upper right", bbox_to_anchor=LEGEND_ANCHOR, frameon=False,
                       handlelength=.7, handletextpad=.3, borderpad=0, labelspacing=.15, borderaxespad=0)
    ticks = [position + offset for position in positions.position_index for offset in INSERT_OFFSETS.values()]
    ax.set_xticks(ticks, [variant.upper() for _ in positions.position_index for variant in INSERT_OFFSETS])
    for position in positions.itertuples():
        for label, value in [("A context", position.context), ("A position", position.position_id[1:])]:
            ax.annotate(value, (position.position_index, 0), xycoords=ax.get_xaxis_transform(),
                        xytext=(0, ANNOTATION_ROWS_PT[label]), textcoords="offset points", ha="center", va="center")
    for label, offset in ANNOTATION_ROWS_PT.items():
        ax.annotate(label, (ROW_LABEL_X, 0), xycoords=ax.get_xaxis_transform(),
                    xytext=(0, offset), textcoords="offset points", ha="right", va="center", annotation_clip=False)
    # Only significant tests are drawn, stacked in points above each site's own data; deepen the axis until the tallest fits.
    rows = rows.assign(top=rows.mean_percent + rows.se_percent)
    significant = tests.loc[tests.stars.ne("ns") & tests.comparison.eq("WT/MUT")].set_index("position_id").stars
    figure.canvas.draw()
    height_pt = ax.get_window_extent().height * 72 / figure.dpi
    level_pt = BRACKET_GAP_PT + BRACKET_LEG_PT + LABEL_HEIGHT_PT
    floor = rows.top.max()
    for _ in range(ITERATIONS):
        ylim = floor
        for _ in range(ITERATIONS):
            per_pt = ylim / height_pt
            brackets = []
            for position in positions.itertuples():
                tops = rows.loc[rows.position_id.eq(position.position_id)].groupby("insert_type").top.max()
                tops = {variant: top for variant, top in tops.items()}
                if position.position_id in significant:
                    brackets.append((position.position_index + INSERT_OFFSETS["wt"], position.position_index + INSERT_OFFSETS["mut"],
                                     max(tops.values()) + per_pt * (level_pt - LABEL_HEIGHT_PT), significant[position.position_id]))
            ylim = max(floor, max([y + per_pt * LABEL_HEIGHT_PT for _, _, y, _ in brackets], default=0))
        ax.set_ylim(0, ylim)
        leg = per_pt * BRACKET_LEG_PT
        added = []
        for left, right, y, stars in brackets:
            added += ax.plot([left, left, right, right], [y - leg, y, y, y - leg], color="black", linewidth=BRACKET_WIDTH, solid_capstyle="butt")
            added.append(ax.annotate(stars, ((left + right) / 2, y), xytext=(0, STAR_OFFSET_PT), textcoords="offset points", ha="center", va="bottom"))
        # Keep the common legend position unless its column would cover a marker, bracket, or star.
        legend.set_bbox_to_anchor(LEGEND_ANCHOR)
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        box = legend.get_window_extent(renderer)
        clearance = LEGEND_CLEARANCE_PT * figure.dpi / 72
        axis_bottom = ax.get_window_extent(renderer).y0
        # Lines hold one point per site, so clear each point rather than the box around the whole series.
        points = [point for line in ax.lines for point in line.get_transform().transform(line.get_xydata())]
        boxes = [text.get_window_extent(renderer) for text in ax.texts]
        data_top = max([y for x, y in points if box.x0 - clearance <= x <= box.x1 + clearance]
                       + [item.y1 for item in boxes if item.y0 >= axis_bottom and item.x1 >= box.x0 - clearance and item.x0 <= box.x1 + clearance],
                       default=box.y0 - clearance)
        shift = max(0, data_top + clearance - box.y0)
        legend.set_bbox_to_anchor((LEGEND_ANCHOR[0], LEGEND_ANCHOR[1] + shift / ax.get_window_extent(renderer).height))
        figure.canvas.draw()
        if legend.get_window_extent(renderer).y1 <= figure.bbox.height:
            break
        # The raised legend leaves the canvas: give the data more room beneath it and redraw the brackets.
        for artist in added:
            artist.remove()
        floor *= HEADROOM_STEP
    require(legend.get_window_extent(renderer).y1 <= figure.bbox.height, f"Distance legend exceeds the canvas in {output.source}")
    print(f"{output.source}: legend raised {shift * 72 / figure.dpi:.2f} pt")
    boxes = [text.get_window_extent() for text in ax.texts]
    require(all(not left.overlaps(right) for i, left in enumerate(boxes) for right in boxes[i + 1:]),
            f"Overlapping distance-figure annotations in {output.source}")
    print(f"{output.source}: paired construct tests")
    print(tests.groupby("comparison").stars.value_counts().unstack(fill_value=0).to_string())
    save_pdf(figure, output.source, png=True)
