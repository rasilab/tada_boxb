"""Recorder sites, with WT/MUT means and SEM across combined or separate placements.

Each construct contributes one editing fraction per site, and a WT entry appearing once per paired
mutant is deduplicated before averaging. SEM describes variation across constructs. Standard recorder
background only, with the four cellular placements retained and no in vitro exclusions imported.
"""

import sys
import os

import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from matplotlib.ticker import MultipleLocator

from analysis_helpers import (annotation, comparison_samples, comparison_set_ids, excluded_constructs,
                              matched_stem_pairs, recorder_positions, require)
from in_vitro_reuse import style as in_vitro_style, helpers as in_vitro_helpers
from plot_style import DISTANCE_SHADES, FONT, FONT_SIZE, MUTANT_COLORS, palettes, plt, save_figure

# Snakemake passes the comparison set so an individual panel can be rebuilt independently.
COMPARISON_SETS = sys.argv[1:] or comparison_set_ids()
BOTTOM_HEIGHT = .4
INSERT_OFFSETS = {"wt": -.16, "mut": .16}
DISTANCE_OFFSETS = [-.075, -.025, .025, .075]
MARKER_SIZE = 2.5
ERROR_WIDTH = .4
ERROR_CAPSIZE = 1
Y_TICK_STEP = 5
CONTEXT_Y = -.20
POSITION_Y = -.38
ROW_LABEL_X = .48
DARK_FACTOR = .7
INTERMEDIATE_WEIGHT = .55
LEGEND_TITLE = "Recorder placement\nrelative to boxB"
Y_LABEL = "% Edited Adenosine"
Y_LABEL_PAD_PT = 3
LABEL_CLEARANCE_PT = .5
LEGEND_LABEL_SPACING = .3
STATISTICS_DIRECTORY = "analysis/boxb_in_vivo_sequencing/data/figure_data"
EXPECTED_SITE_TESTS = 24
BRACKET_GAP_PT = 3
BRACKET_LEG_PT = 2
LABEL_HEIGHT_PT = 9
STAR_OFFSET_PT = -2
BRACKET_WIDTH = .4
ITERATIONS = 20

require(set(COMPARISON_SETS) <= set(comparison_set_ids()), "Unknown comparison set")
panel_selections = annotation("recorder_site_panels").set_index("comparison_set", verify_integrity=True)
require(set(COMPARISON_SETS) <= set(panel_selections.index), "Missing recorder-site panel selection")
require(excluded_constructs().empty, "Flagged placements require an explicit with/without comparison")
labels, light, colors = palettes()
positions = recorder_positions().sort_values("sense_position").reset_index(drop=True)
positions["position_index"] = np.arange(1, len(positions) + 1)
measurements = positions.assign(sites=positions.position_id, label=positions.position_id,
                                context=positions.derived_context)[
    ["position_id", "sites", "label", "context", "editing_column"]]

for set_id in COMPARISON_SETS:
    panel = panel_selections.loc[set_id]
    selected = comparison_samples(set_id)
    sample_ids = panel.sample_ids.split(";")
    require(len(sample_ids) == len(set(sample_ids)) and set(sample_ids) <= set(selected.sample_id),
            "Repeated or unannotated recorder-site library")
    selected = selected.set_index("sample_id").loc[sample_ids].reset_index()
    require(len(selected) > 0, "The selected comparison has no libraries")
    matched = matched_stem_pairs(selected, measurements)
    matched["percent_edited"] = 100 * matched.edited_umis / matched.callable_umis
    statistics = in_vitro_helpers.recorder_site_tests(matched)
    require(statistics.groupby("sample_id").size().eq(EXPECTED_SITE_TESTS).all(), "Incomplete site tests")
    os.makedirs(STATISTICS_DIRECTORY, exist_ok=True)
    statistics.to_csv(f"{STATISTICS_DIRECTORY}/site_placement_{set_id}_statistics.csv", index=False)
    pairs = matched.drop_duplicates(["sample_id", "position_id", "barcode", "insert"])
    require(pairs.callable_umis.gt(0).all(), "A construct and recorder site carries no callable molecule")
    geometry_keys = ["target_pos_to_boxb", "target_dist"] if panel.split_distances else []
    summary = pairs.groupby(["sample_id", "enzyme", "panel_order", "position_id", *geometry_keys,
                             "insert_type"], as_index=False).agg(
        mean_percent=("percent_edited", "mean"), se_percent=("percent_edited", "sem"),
        constructs=("insert", "size"))
    require(summary.se_percent.notna().all(), "A plotted mean has insufficient constructs for SEM")
    summary = summary.merge(positions[["position_id", "position_index", "derived_context"]],
                            on="position_id", validate="many_to_one")
    geometries = pairs[["target_pos_to_boxb", "target_dist"]].drop_duplicates().sort_values(
        "target_dist").reset_index(drop=True)
    require(len(geometries) == 4, "Expected four designed recorder placements")
    sites = positions.position_id.tolist()
    require(len(summary) == len(selected) * len(sites) * (len(geometries) if panel.split_distances else 1) * len(INSERT_OFFSETS),
            "Incomplete recorder site/placement summary")
    if panel.split_distances:
        in_vitro_labels, _, _ = in_vitro_style.palettes()
        distance_colors = {in_vitro_labels[enzyme]: shades
                           for enzyme, shades in in_vitro_style.distance_palettes().items()}

    figure, axes = plt.subplots(len(selected), 1,
                                figsize=(panel.width, panel.panel_height * len(selected) + BOTTOM_HEIGHT),
                                layout="constrained", sharex=True, sharey=True, squeeze=False)
    y_max = np.ceil((summary.mean_percent + summary.se_percent).max() / Y_TICK_STEP) * Y_TICK_STEP
    for row, library in enumerate(selected.itertuples()):
        ax = axes[row][0]
        if panel.split_distances:
            base, pale = np.array(to_rgb(colors[library.enzyme])), np.array(to_rgb(light[library.enzyme]))
            wt_colors = distance_colors.get(labels[library.enzyme],
                                            [base * DARK_FACTOR, base,
                                             base * (1 - INTERMEDIATE_WEIGHT) + pale * INTERMEDIATE_WEIGHT, pale])
        for insert_type, offset in INSERT_OFFSETS.items():
            data = summary.loc[summary.sample_id.eq(library.sample_id) & summary.insert_type.eq(insert_type)]
            if panel.split_distances:
                for index, geometry in enumerate(geometries.itertuples()):
                    points = data.loc[data.target_pos_to_boxb.eq(geometry.target_pos_to_boxb)
                                      & data.target_dist.eq(geometry.target_dist)].sort_values("position_index")
                    ax.errorbar(points.position_index + offset + DISTANCE_OFFSETS[index], points.mean_percent,
                                yerr=points.se_percent, fmt="o", markersize=MARKER_SIZE, linewidth=ERROR_WIDTH,
                                capsize=ERROR_CAPSIZE,
                                color=wt_colors[index] if insert_type == "wt" else DISTANCE_SHADES[index])
            else:
                points = data.sort_values("position_index")
                ax.errorbar(points.position_index + offset, points.mean_percent, yerr=points.se_percent,
                            fmt="o", markersize=MARKER_SIZE, linewidth=ERROR_WIDTH, capsize=ERROR_CAPSIZE,
                            color=colors[library.enzyme] if insert_type == "wt" else MUTANT_COLORS["dark"],
                            label=insert_type.upper())
        ax.set(ylim=(0, y_max), xlim=(.55, len(sites) + .45))
        ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_STEP))
        ax.set_xticks(positions.position_index, [])
        ax.tick_params(axis="x", length=0)

    # Draw only the pre-specified WT/MUT contrast. The placement comparisons stay in the
    # statistics CSV because this panel does not separate placements along its axis.
    figure.canvas.draw()
    height_pt = axes[0][0].get_window_extent().height * 72 / figure.dpi
    level_pt = BRACKET_GAP_PT + BRACKET_LEG_PT + LABEL_HEIGHT_PT
    tops = summary.assign(top=summary.mean_percent + summary.se_percent)
    marked = {}
    for row, library in enumerate(selected.itertuples()):
        at_site = statistics.loc[statistics.sample_id.eq(library.sample_id)
                                 & statistics.comparison.eq("WT/MUT")]
        require(len(at_site) == len(positions), "Every recorder site needs one WT/MUT test")
        marked[row] = at_site.loc[at_site.stars.ne("ns")].set_index("position_id").stars
        print(f"{library.sample_id}: WT/MUT significant at {len(marked[row])} of {len(at_site)} sites")

    def bracket_tops(row, library, per_pt):
        """Each drawn bracket sits above its own site's tallest error bar."""
        rows_here = tops.loc[tops.sample_id.eq(library.sample_id)]
        placed = []
        for position in positions.itertuples():
            if position.position_id not in marked[row].index:
                continue
            at = rows_here.loc[rows_here.position_id.eq(position.position_id)].top.max()
            placed.append((position, at + per_pt * (level_pt - LABEL_HEIGHT_PT)))
        return placed

    # The bracket height depends on the axis limit, which the brackets in turn set.
    ylim = y_max
    for _ in range(ITERATIONS):
        per_pt = ylim / height_pt
        needed = [y_max]
        for row, library in enumerate(selected.itertuples()):
            needed += [y + per_pt * LABEL_HEIGHT_PT for _, y in bracket_tops(row, library, per_pt)]
        ylim = max(needed)
    ylim = np.ceil(ylim / Y_TICK_STEP) * Y_TICK_STEP
    per_pt, leg = ylim / height_pt, (ylim / height_pt) * BRACKET_LEG_PT
    for row, library in enumerate(selected.itertuples()):
        ax = axes[row][0]
        ax.set_ylim(0, ylim)
        for position, y in bracket_tops(row, library, per_pt):
            left = position.position_index + INSERT_OFFSETS["wt"]
            right = position.position_index + INSERT_OFFSETS["mut"]
            ax.plot([left, left, right, right], [y - leg, y, y, y - leg], color="black",
                    linewidth=BRACKET_WIDTH, solid_capstyle="butt")
            ax.annotate(marked[row][position.position_id], ((left + right) / 2, y),
                        xytext=(0, STAR_OFFSET_PT), textcoords="offset points", ha="center", va="bottom")

    ax = axes[-1][0]
    for position in positions.itertuples():
        context = HPacker(children=[TextArea(base, textprops=dict(fontfamily=FONT, fontsize=FONT_SIZE,
                           fontweight="bold" if index == 1 else "normal"))
                           for index, base in enumerate(position.derived_context)], align="baseline", pad=0, sep=0)
        ax.add_artist(AnnotationBbox(context, (position.position_index, CONTEXT_Y),
                                    xycoords=ax.get_xaxis_transform(), box_alignment=(.5, .5),
                                    frameon=False, pad=0, annotation_clip=False))
        ax.text(position.position_index, POSITION_Y, str(position.sense_position), ha="center",
                va="center", transform=ax.get_xaxis_transform())
    for y, label in [(CONTEXT_Y, "A context"), (POSITION_Y, "A position")]:
        ax.text(ROW_LABEL_X, y, label, ha="right", va="center", transform=ax.get_xaxis_transform())
    axes[0, 0].set_ylabel(Y_LABEL, labelpad=Y_LABEL_PAD_PT)
    if panel.split_distances:
        handles = [Line2D([], [], marker="o", linestyle="none", color=DISTANCE_SHADES[index],
                          markersize=MARKER_SIZE, label=f"{geometry.target_pos_to_boxb}′, {geometry.target_dist} nt")
                   for index, geometry in enumerate(geometries.itertuples())]
        figure.legend(handles=handles, loc="outside right center", frameon=False, title=LEGEND_TITLE,
                      handlelength=.7, borderpad=0, labelspacing=.4)
    else:
        handles = [Line2D([], [], marker="o", linestyle="none", markersize=MARKER_SIZE,
                          color=colors[selected.iloc[0].enzyme] if variant == "wt" else MUTANT_COLORS["dark"],
                          label=variant.upper()) for variant in INSERT_OFFSETS]
        axes[0, 0].legend(handles=handles, loc="upper right", frameon=False,
                         handlelength=.7, borderpad=0, labelspacing=LEGEND_LABEL_SPACING)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    for ax in axes.flat:
        boxes = [text.get_window_extent(renderer) for text in ax.get_xticklabels()]
        require(all(right.x0 - left.x1 >= LABEL_CLEARANCE_PT * figure.dpi / 72
                    for left, right in zip(boxes, boxes[1:])), "WT/MUT labels overlap")
    save_figure(figure, f"site_placement_{set_id}")

    print(f"\n{set_id} batch, editing by recorder site, split_distances={panel.split_distances}, "
          "construct means (%), error bars are SEM:")
    wide = summary.pivot(index=["sample_id", *geometry_keys, "insert_type"],
                         columns="position_id", values="mean_percent")[sites]
    print(wide.to_string(float_format="%.2f"))
    print(f"  each mean uses {int(summary.constructs.min())}-{int(summary.constructs.max())} distinct constructs; "
          "SEM does not represent biological replication")
    print("  no in vivo placements are excluded; with/without-exclusion estimates are identical")
