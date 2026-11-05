"""Pooled matched-WT effects; sites share molecules and are not reaction replicates."""

import math
import os
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import ScaledTranslation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import gaussian_kde
from analysis_helpers import annotation, MEASUREMENT_SEPARATOR, designed_boxb_sequences, designed_variant_effects, drawn_at, require, single_mutant_data, substitution_measurements, boxb_wildtype, reaction_set, reaction_set_ids, scale_groups, RNA_BASES
from plot_style import plt, palettes, save_pdf, logo_letter, place_shared_ylabel, centre_shared_xlabel, LABEL_GAP

# One page geometry per scale group. The editors set draws its heat map cells, and the logo glyphs
# under them, a tenth narrower than the concentration series does, so its block fits a narrower column
# of an assembled figure; the margins stay the same width in inches, which keeps the central row names
# level with the logo label above them and leaves room for the measurement named at the right.
# Each group shares a base logo scale; per-reaction presentation overrides live in the CSV.
GEOMETRY = {"lambdaN_tada": dict(width=5.1, left=.153, right=.941, row_label_x=-.9, ylabel_x=-.072,
                                 logo_height=.64, panel_gap=.10),
            "editors": dict(width=4.7, left=.166, right=.936, row_label_x=-1.0, ylabel_x=-.080,
                            logo_height=.45, panel_gap=.10),
            "time": dict(width=5.1, left=.153, right=.941, row_label_x=-.9, ylabel_x=-.072,
                         logo_height=.64, panel_gap=.10)}
TILE_WIDTH_SCALE = .925
# Narrow tiles and aligned logos while keeping margins and label offsets fixed in inches.
for geometry in GEOMETRY.values():
    panel_width = geometry["width"] * (geometry["right"] - geometry["left"])
    figure_width = geometry["width"] - panel_width * (1 - TILE_WIDTH_SCALE)
    geometry["left"] *= geometry["width"] / figure_width
    geometry["right"] = 1 - (1 - geometry["right"]) * geometry["width"] / figure_width
    geometry["width"] = figure_width
    geometry["ylabel_x"] /= TILE_WIDTH_SCALE
    geometry["row_label_x"] = -.5 + (geometry["row_label_x"] + .5) / TILE_WIDTH_SCALE
PANEL_HEIGHT = .64
CENTER_HEIGHT = .73
TOP_MARGIN = .06
BOTTOM_MARGIN = .16
TITLE_HEIGHT = .18
RIGHT_LABEL_GAP = .04
# Narrowing the figure around two square panels is what brings them together: the equal aspect fixes
# the panel size, so the width the panels do not use is white space between them. Both sets take one
# panel size; compact output removes the editor label's released left margin.
CORRELATION_LAYOUTS = {
    "lambdaN_tada": dict(size=(2.36, 1.31), ylabel_gap=LABEL_GAP, compact=False, png=False),
    "editors": dict(size=(2.36, 1.31), ylabel_gap=1 / 72, compact=True, png=True),
    "time": dict(size=(2.36, 1.31), ylabel_gap=LABEL_GAP, compact=False, png=True),
}
# Label every second power of two where they fit, and thin the lattice where they do not: an axis deep
# enough to hold the weakest designed variants carries its labels closer together than a panel this
# size can separate them. The step is chosen against the drawn numbers rather than assumed, since the
# panel size, the axis depth and the label widths together decide it.
# A comparison may be stacked for the reaction sets its annotation names: two square panels one
# above the other rather than side by side, which is a third of the width and fits a figure that
# sets three panels across. The same comparison stays side by side wherever it is not named.
STACKED_SIZE = (1.52, 2.62)
CORRELATION_TICK_STEPS = [2, 3, 4]
# Numbers that merely fail to overlap still read as one number at this size: 1.56 and 6.25 a scant
# point apart read as 1.566.25. A step is accepted only once its numbers stand this far apart.
TICK_LABEL_GAP = 2
LOWEST_TICK = -4
# Each comparison names the variant set it draws: the 57 single substitutions, or every sequence the
# combinatorial library designs, singles included. The single-substitution panels may leave one point
# below the axis rather than deepen both panels for it; the designed panels carry hundreds of points,
# where one deep point costs the whole set nothing, so their axis reaches every drawn variant.
# Those panels also draw an order of magnitude more points into the same square, so the markers are
# smaller and translucent, which keeps the dense centre of the distribution readable. They are read on
# their own rather than assembled over another panel, so they carry a white page instead of a
# transparent one, on which translucent markers cannot be read against whatever lies behind them.
VARIANT_POINTS_OUTSIDE = {"single": 1, "designed": 0}
VARIANT_POINT_STYLES = {"single": dict(s=7, alpha=1), "designed": dict(s=1.6, alpha=.45)}
VARIANT_TRANSPARENT = {"single": True, "designed": False}
STATISTICS_PATH = "analysis/boxb_in_vitro_sequencing/data/figure_data/boxb_substitution_correlation_statistics.csv"
COLOR_STEP = .5
WT_COLOR = "#E8E8E8"
REFERENCE_COLOR = "#999999"
ANNOTATION_COLOR = "#000000"
NUMBER_ON_LIGHT = "#6B6B6B"
NUMBER_ON_DARK = "#EDEDED"
NUMBER_SHADE_THRESHOLD = .6
SMALL_SIZE = 6
BASE_OFFSET = 16
POSITION_OFFSET = 40
STRUCTURE_OFFSET = 25.5
STRUCTURE_BAR = 4
STRUCTURE_GAP = .25
DENSITY_POINTS = 256
DENSITY_WIDTH_RATIO = .16
DENSITY_GAP_IN = .04
DENSITY_ALPHA = .6
LOGO_AXIS_STEP = 10
LOGO_TICKS = [0, -50]
LOGO_LETTER_WIDTH = .78
LOGO_YLABEL = "Editing loss\n% of wild type"
ZERO_LINE_WIDTH = .4
TOLERANCE = 1e-9

all_measurements = substitution_measurements()
logo_geometry = annotation("boxb_reaction_sets").set_index(["reaction_set", "sample_id"])
logo_axis_minima = logo_geometry.logo_axis_min
logo_magnifications = logo_geometry.logo_magnification.fillna(1)
require(logo_axis_minima.dropna().le(min(LOGO_TICKS)).all() and np.isfinite(logo_axis_minima.dropna()).all(), "A declared logo axis must contain both labelled ticks")
require(logo_magnifications.gt(0).all() and np.isfinite(logo_magnifications).all(), "Invalid logo magnification")
all_correlations = annotation("boxb_substitution_correlations")
# A comparison names the sets it stacks, or none of them.
all_correlations["stacked_sets"] = all_correlations.stacked_sets.fillna("") if "stacked_sets" in all_correlations else ""
require(len(all_correlations) > 0 and all_correlations.notna().all().all() and not all_correlations.source.duplicated().any() and not all_correlations.duplicated(["x_position_id", "y_position_id", "variant_set"]).any(), "Invalid substitution correlation selections")
require(all_correlations[["x_position_id", "y_position_id"]].isin(all_measurements.position_id.tolist()).all().all() and all_correlations.x_position_id.ne(all_correlations.y_position_id).all(), "Invalid correlation measurements")
require(all_correlations.variant_set.isin(VARIANT_POINTS_OUTSIDE).all(), "Every comparison must name a known variant set")
groups = scale_groups()
group_of = {set_id: group for group, members in groups.items() for set_id in members}
set_ids = reaction_set_ids()
labels, _, colors = palettes()
regions, wildtype, _ = boxb_wildtype()
position_count = len(wildtype)


def shifted_below(figure, ax, points):
    """Place annotations a fixed point distance below one axis, independent of panel height."""
    return ax.get_xaxis_transform() + ScaledTranslation(0, -points / 72, figure.dpi_scale_trans)


def annotate_boxb_positions(figure, ax, row_label_x):
    """Label wild-type bases, stem-loop regions, and boxB positions once, between the two enzyme blocks."""
    for x, base in enumerate(wildtype):
        ax.text(x, 0, base, ha="center", fontsize=SMALL_SIZE, color=ANNOTATION_COLOR, transform=shifted_below(figure, ax, BASE_OFFSET))
        ax.text(x, 0, str(x + 1), ha="center", color=ANNOTATION_COLOR, transform=shifted_below(figure, ax, POSITION_OFFSET))
    # Both row names break onto two lines, which keeps the left margin narrow enough that naming the
    # measurement beside each heat map does not widen the page.
    for points, label in [(BASE_OFFSET, "Wild-type\nbase"), (POSITION_OFFSET, "boxB\nposition")]:
        ax.text(row_label_x, 0, label, ha="right", va="center", multialignment="right", color=ANNOTATION_COLOR, transform=shifted_below(figure, ax, points))
    # Each region label sits on its own span line, which breaks around the measured width of the text.
    row = shifted_below(figure, ax, STRUCTURE_OFFSET)
    renderer = figure.canvas.get_renderer()
    for region in regions.itertuples():
        left, right = region.start - 1.4, region.end - .6
        label = ax.text((left + right) / 2, 0, region.label, ha="center", va="center", color=ANNOTATION_COLOR, transform=row)
        box = label.get_window_extent(renderer)
        start, end = (ax.transData.inverted().transform((x, 0))[0] for x in [box.x0, box.x1])
        require(left < start - STRUCTURE_GAP and end + STRUCTURE_GAP < right, f"Region label does not fit its span: {region.label}")
        for span in [[left, start - STRUCTURE_GAP], [end + STRUCTURE_GAP, right]]:
            ax.plot(span, [0, 0], color=ANNOTATION_COLOR, linewidth=.4, transform=row, clip_on=False)
        ax.plot([left, right], [0, 0], linestyle="none", marker="|", markersize=STRUCTURE_BAR, markeredgewidth=.4, color=ANNOTATION_COLOR, transform=row, clip_on=False)


def percent_ticks(limits, step):
    """Powers-of-two ticks on a lattice anchored on wild type, so 100% always carries a label.

    The lowest label stays at 6.25% of wild type wherever the axis itself ends, so a panel on an axis
    deepened to hold the weakest substitutions carries the same labels as the panels it is read against.
    An axis that reaches deeper than that is labelled down to its own end rather than leaving its lower
    half blank. Anchoring on wild type rather than on the lowest label is what lets the step change
    without moving the labels the panels are read by: 100% is a tick at every step.
    """
    lower, upper = limits
    return list(range(step * math.ceil(min(LOWEST_TICK, lower) / step), step * math.floor(upper / step) + 1, step))


def logo_stacks(selected, enzyme_count):
    """Mean loss of the three substitutions relative to matched WT, drawn below zero."""
    activities = selected.pivot(index=["sample_id", "boxb_position"], columns="nucleotide", values="pooled_ratio_to_matched_wt")[list(RNA_BASES)]
    wildtype_mask = selected.pivot(index=["sample_id", "boxb_position"], columns="nucleotide", values="is_wildtype")[list(RNA_BASES)]
    require(activities.shape == (enzyme_count * position_count, len(RNA_BASES)) and np.isfinite(activities.to_numpy()).all() and activities.gt(0).all().all(), "Incomplete substitution activity matrix")
    require(wildtype_mask.sum(axis=1).eq(1).all() and activities.where(wildtype_mask).max(axis=1).eq(1).all(), "Wild-type reference is not one activity of one per position")
    # Restore the original three-mutant mean and WT reference, with loss below zero.
    # Positions whose average activity equals or exceeds WT have no downward loss stack;
    # their unmodified per-substitution activities remain visible in the heat map.
    mutant_means = activities.mask(wildtype_mask).mean(axis=1)
    changes = (100 * (mutant_means - 1)).clip(upper=0)
    weights = activities.div(activities.sum(axis=1), axis=0)
    heights = weights.mul(-changes, axis=0)
    require(changes.between(-100, 0).all() and activities.mask(wildtype_mask).notna().sum(axis=1).eq(len(RNA_BASES) - 1).all(), "Invalid matched-WT editing loss")
    require(np.allclose(weights.sum(axis=1), 1) and np.allclose(heights.sum(axis=1), -changes), "Invalid logo weights or letter heights")
    return heights, changes


# Build every set before drawing any, so the sets of one scale group share a scale per measurement
# and per comparison. A colour, a stack height, and an axis then mean the same thing across the sets
# that will sit beside each other, and the shared scale is the widest any of them needs, so nothing
# is clipped. Sets of different groups are scaled apart: a figure built from one should not be
# stretched by a set it will never sit beside.
datasets = {set_id: single_mutant_data(set_id)[0] for set_id in set_ids}
drawn_measurements = {set_id: drawn_at(all_measurements, set_id) for set_id in set_ids}
drawn_correlations = {set_id: drawn_at(all_correlations, set_id) for set_id in set_ids}
for set_id, tiles in datasets.items():
    require(tiles.wildtype_rna.eq(wildtype).all(), "Tiles carry a different wild-type sequence")
    # A set draws only the measurements and comparisons its own rows declare, so a measurement whose
    # rarest sequences reach zero edited molecules in that set is left out by annotation.
    require(drawn_correlations[set_id][["x_position_id", "y_position_id"]].isin(drawn_measurements[set_id].position_id.tolist()).all().all(), f"A comparison in {set_id} needs a measurement that set does not draw")

# The whole designed library is read only for the sets whose comparisons ask for it, and only at the
# recorder sites those comparisons name, since each site is a second pass over every reaction's counts.
designed_count = len(designed_boxb_sequences())
designed_datasets = {}
for set_id in set_ids:
    chosen = drawn_correlations[set_id]
    chosen = chosen.loc[chosen.variant_set.eq("designed")]
    if len(chosen):
        designed_datasets[set_id] = designed_variant_effects(set_id, sorted(set(chosen.x_position_id) | set(chosen.y_position_id)))


def comparison_points(set_id, comparison, reaction):
    """One enzyme's paired effects for a comparison, as the variant set it names defines them.

    Wild type is left out of both sets: it is the reference every ratio is taken against, so it would
    enter as one point at 100% on both axes by construction. A designed sequence with no edited
    molecule at either site has no logarithm to place on these axes; it is dropped here and counted,
    rather than given a floor value that would read as a measurement.
    """
    sites = [comparison.x_position_id, comparison.y_position_id]
    if comparison.variant_set == "single":
        tiles = datasets[set_id]
        selected = tiles.loc[tiles.position_id.isin(sites) & ~tiles.is_wildtype & tiles.sample_id.eq(reaction.sample_id)]
        points = selected.pivot(index=["boxb_position", "nucleotide"], columns="position_id", values="log2_pooled_ratio_to_matched_wt")
        require(len(points) == (len(RNA_BASES) - 1) * position_count and np.isfinite(points.to_numpy()).all(), "Incomplete paired substitutions")
        return points[sites], 0
    variants = designed_datasets[set_id]
    selected = variants.loc[variants.position_id.isin(sites) & variants.n_substitutions.gt(0) & variants.sample_id.eq(reaction.sample_id)]
    points = selected.pivot(index="hairpin", columns="position_id", values="log2_pooled_ratio_to_matched_wt")
    require(len(points) == designed_count - 1, "Incomplete paired designed variants")
    finite = np.isfinite(points).all(axis=1)
    return points.loc[finite, sites], int((~finite).sum())

color_limits, logo_limits, scatter_limits = {}, {}, {}
for group, members in groups.items():
    for measurement in all_measurements.itertuples():
        drawing = [set_id for set_id in members if measurement.position_id in drawn_measurements[set_id].position_id.tolist()]
        if not drawing:
            continue
        selections = [datasets[set_id].loc[datasets[set_id].position_id.eq(measurement.position_id)] for set_id in drawing]
        color_limits[(group, measurement.position_id)] = max([COLOR_STEP] + [math.ceil(-selected.log2_pooled_ratio_to_matched_wt.min() / COLOR_STEP) * COLOR_STEP for selected in selections])
        if measurement.include_logo:
            peak = max(-logo_stacks(selected, selected.sample_id.nunique())[1].min() for selected in selections)
            logo_limits[group] = max(logo_limits.get(group, -min(LOGO_TICKS)), math.ceil(peak / LOGO_AXIS_STEP) * LOGO_AXIS_STEP)
    for comparison in all_correlations.itertuples():
        values = []
        for set_id in members:
            if comparison.source not in drawn_correlations[set_id].source.tolist():
                continue
            for reaction in reaction_set(set_id).itertuples():
                values.append(comparison_points(set_id, comparison, reaction)[0].to_numpy().ravel())
        if not values:
            continue
        combined = np.concatenate(values)
        require(np.isfinite(combined).all(), f"A substitution has no ratio to place on a log axis in {comparison.source}")
        # The axis reaches every substitution but the one furthest below the rest, which would stretch
        # both panels of its group on its own; POINTS_OUTSIDE fixes how many may sit outside, and the
        # rest of the distribution is drawn rather than pressed against a floor. The labelled ticks stay
        # where they are, so a panel on a deeper axis still reads against the others.
        ordered = np.sort(combined)
        allowed = VARIANT_POINTS_OUTSIDE[comparison.variant_set]
        lower = min(LOWEST_TICK, math.floor(ordered[allowed]))
        outside = ordered[ordered < lower]
        require(len(outside) <= allowed, f"More than {allowed} substitution(s) would sit below the axis in {comparison.source}")
        if len(outside):
            print(f"{group} {comparison.source}: {len(outside)} substitution(s) at "
                  + ", ".join(f"{100 * 2. ** value:.2f}%" for value in outside)
                  + f" of wild type sit below the {100 * 2. ** lower:.2f}% axis end; every other substitution is drawn")
        scatter_limits[(group, comparison.source)] = (lower, math.ceil(combined.max()))
for group in groups:
    print(f"Shared scales for {group}: "
          + ", ".join(f"{position_id} colour +/-{limit:.2f}" for (name, position_id), limit in color_limits.items() if name == group)
          + f", logo scale {logo_limits[group] / GEOMETRY[group]['logo_height']:g} percentage points of WT per inch"
          + ", " + ", ".join(f"{source} axis {100 * 2. ** lower:g}-{100 * 2. ** upper:g}%" for (name, source), (lower, upper) in scatter_limits.items() if name == group))

# Every figure of this family is drawn once per reaction set, which only its filename records.
correlation_statistics = []
for set_id in set_ids:
    group = group_of[set_id]
    metadata = reaction_set(set_id)
    # Panels are keyed by their reaction, since a set may draw one enzyme twice; the title tells
    # them apart and the colour still follows the enzyme.
    titles = {row.sample_id: row.panel_title if isinstance(row.panel_title, str) else labels[row.tada_type] for row in metadata.itertuples()}
    enzyme_of = dict(zip(metadata.sample_id, metadata.tada_type))
    geometry = GEOMETRY[group]
    tiles = datasets[set_id]
    measurements, correlations = drawn_measurements[set_id], drawn_correlations[set_id]
    # Each enzyme gets its heat map with its logo directly below; both blocks share the central position labels.
    for measurement in measurements.loc[measurements.heatmap_source.notna()].itertuples():
        selected = tiles.loc[tiles.position_id.eq(measurement.position_id)]
        limit = color_limits[(group, measurement.position_id)]
        panel_kinds = ["heatmap", "logo"] if measurement.include_logo else ["heatmap"]
        if "logo" in panel_kinds:
            heights, changes = logo_stacks(selected, len(metadata))
            logo_limit = logo_limits[group]
        # Shorter axes retain the base percentage per inch unless magnification is declared.
        # Zero stays below the heat map; each label follows its panel's centre.
        panel_heights, panel_limits, block_heights = {}, {}, {}
        for reaction in metadata.itertuples():
            enzyme = reaction.sample_id
            panel_heights[enzyme] = {"heatmap": PANEL_HEIGHT}
            if "logo" in panel_kinds:
                declared_minimum = logo_axis_minima.loc[(set_id, reaction.sample_id)]
                panel_limits[enzyme] = -declared_minimum if np.isfinite(declared_minimum) else logo_limit
                panel_heights[enzyme]["logo"] = geometry["logo_height"] * panel_limits[enzyme] / logo_limit * logo_magnifications.loc[(set_id, reaction.sample_id)]
            block_heights[enzyme] = sum(panel_heights[enzyme].values()) + (len(panel_kinds) - 1) * geometry["panel_gap"]
        # Each block carries its enzyme's name above its heat map, so the block needs a title band of
        # its own and the page needs no heading over the whole figure.
        figure_height = TOP_MARGIN + 2 * TITLE_HEIGHT + sum(block_heights.values()) + CENTER_HEIGHT + BOTTOM_MARGIN
        figure = plt.figure(figsize=(geometry["width"], figure_height))
        panels, shared = {}, None
        for block, enzyme in enumerate(metadata.sample_id):
            block_top = TOP_MARGIN + TITLE_HEIGHT + block * (TITLE_HEIGHT + CENTER_HEIGHT) + sum(block_heights[previous] for previous in metadata.sample_id.iloc[:block])
            for index, kind in enumerate(panel_kinds):
                height = panel_heights[enzyme][kind]
                top = block_top + sum(panel_heights[enzyme][above] + geometry["panel_gap"] for above in panel_kinds[:index])
                ax = figure.add_axes([geometry["left"], 1 - (top + height) / figure_height, geometry["right"] - geometry["left"], height / figure_height], sharex=shared)
                shared = shared or ax
                panels[(enzyme, kind)] = ax
        for (enzyme, kind), ax in panels.items():
            if kind == "heatmap":
                rows = selected.loc[selected.sample_id.eq(enzyme)]
                values = rows.pivot(index="nucleotide", columns="boxb_position", values="log2_pooled_ratio_to_matched_wt").loc[list(RNA_BASES)].to_numpy()
                mask = rows.pivot(index="nucleotide", columns="boxb_position", values="is_wildtype").loc[list(RNA_BASES)].to_numpy()
                # Colour on the log₂ ratio, which spaces the gradient evenly, but label with the percentage a reader can act on.
                percent = 100 * rows.pivot(index="nucleotide", columns="boxb_position", values="pooled_ratio_to_matched_wt").loc[list(RNA_BASES)].to_numpy()
                require(np.allclose(percent[mask], 100) and np.isfinite(percent).all(), "Wild-type tiles must read 100 percent")
                cmap = LinearSegmentedColormap.from_list("enzyme", [colors[enzyme_of[enzyme]], "white"])
                cmap.set_bad(WT_COLOR)
                ax.set_facecolor(WT_COLOR)
                ax.pcolormesh(np.arange(position_count + 1) - .5, np.arange(5) - .5, np.ma.array(values, mask=mask), cmap=cmap, vmin=-limit, vmax=0, rasterized=False)
                # Recede the numbers against their tile so they stay readable up close without dominating.
                for (y, x), value in np.ndenumerate(values):
                    shade = 0 if mask[y, x] else value / -limit
                    ax.text(x, y, f"{percent[y, x]:.0f}", ha="center", va="center", fontsize=SMALL_SIZE, color=NUMBER_ON_DARK if shade > NUMBER_SHADE_THRESHOLD else NUMBER_ON_LIGHT)
                ax.set(yticks=range(4), yticklabels=list(RNA_BASES), ylabel="Mutant base", ylim=(3.5, -.5), title=titles[enzyme])
                ax.tick_params(length=0, labelsize=SMALL_SIZE)
                # The enzyme names its block from above and the measurement stands beside the heat map it
                # labels, repeated once per block, so each block reads on its own.
                box = ax.get_position()
                figure.text(geometry["right"] + RIGHT_LABEL_GAP / geometry["width"], box.y0 + box.height / 2,
                            f"{measurement.axis_label},\n% of wild type", rotation=270, va="center", ha="left", multialignment="center")
                for spine in ax.spines.values():
                    spine.set_visible(False)
            else:
                for position, stack in heights.loc[enzyme].iterrows():
                    bottom = 0
                    # Upright letters descend from zero, with the largest nearest the baseline.
                    for base, height in stack.sort_values(ascending=False, kind="stable").items():
                        bottom -= height
                        if height <= TOLERANCE:
                            continue
                        logo_letter(ax, base, position - 1, bottom, LOGO_LETTER_WIDTH, height)
                    require(np.isclose(bottom, changes.loc[(enzyme, position)]), "A drawn logo does not match its editing loss")
                scores = changes.loc[enzyme]
                require(scores.min() >= -panel_limits[enzyme] - TOLERANCE, f"The logo axis for {titles[enzyme]} cuts off a stack in {set_id}")
                # Two lines centred on each other, which keeps the label's own column narrow: the second
                # line grows leftward from the shared vertical line, not past the central row names.
                ax.set(ylabel=LOGO_YLABEL, ylim=(-panel_limits[enzyme], 0), yticks=LOGO_TICKS, yticklabels=[f"{value:g}" for value in LOGO_TICKS])
                ax.axhline(0, color=REFERENCE_COLOR, linewidth=ZERO_LINE_WIDTH, clip_on=False)
                ax.spines["bottom"].set_visible(False)
                ax.yaxis.label.set_multialignment("center")
                ax.tick_params(bottom=False)
                print(f"{set_id} {measurement.position_id} logo, {titles[enzyme]}: minimum={scores.min():.10f}% of WT at position {scores.idxmin()}, median={scores.median():.10f}, no-loss positions={sorted(scores.index[scores.eq(0)])}, axis=-{panel_limits[enzyme]:g} to 0, scale={logo_limit / geometry['logo_height'] / logo_magnifications.loc[(set_id, enzyme)]:.2f} percentage points per inch")
            ax.set_xlim(-.5, position_count - .5)
            ax.tick_params(labelbottom=False)
            # Fix every axis label at one offset so the heat map and logo labels share a vertical line.
            ax.yaxis.set_label_coords(geometry["ylabel_x"], .5)
        first_enzyme = metadata.sample_id.iloc[0]
        annotate_boxb_positions(figure, panels[(first_enzyme, panel_kinds[-1])], geometry["row_label_x"])
        if "logo" in panel_kinds:
            units_per_inch = [np.ptp(ax.get_ylim()) / (ax.get_position().height * figure_height) for (enzyme, kind), ax in panels.items() if kind == "logo"]
            expected_scales = [logo_limit / geometry["logo_height"] / logo_magnifications.loc[(set_id, sample)] for sample in metadata.sample_id]
            require(np.allclose(units_per_inch, expected_scales), "Logo scales must match their declared magnifications")
            gaps = [(panels[(sample, "heatmap")].get_position().y0 - panels[(sample, "logo")].get_position().y1) * figure_height for sample in metadata.sample_id]
            require(np.allclose(gaps, geometry["panel_gap"]), "Every logo zero must stay the same short distance below its heat map")
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        for ax in panels.values():
            label_box = ax.yaxis.label.get_window_extent(renderer)
            require(label_box.x0 >= 0 and label_box.y0 >= 0 and label_box.y1 <= figure.bbox.height, "An axis label is outside the figure")
            require(not any(label_box.overlaps(tick.get_window_extent(renderer)) for tick in ax.get_yticklabels()), "An axis label overlaps its tick numbers")
            require(not any(label_box.overlaps(text.get_window_extent(renderer)) for other in panels.values() for text in other.texts), "An axis label overlaps a panel annotation")
        save_pdf(figure, f"{measurement.heatmap_source}_{set_id}", transparent=True, png=True)

    # Plot the heat map value, mutant editing as a percentage of wild type, on axes shared by both enzymes.
    axis_labels = {row.position_id: f"{row.axis_label}, % of wild type" for row in measurements.itertuples()}
    # The y label sits in the narrow strip beside the panels, so it breaks at its own comma.
    stacked_labels = {row.position_id: f"{row.axis_label},\n% of wild type" for row in measurements.itertuples()}
    for comparison in correlations.itertuples():
        axis_limits = scatter_limits[(group, comparison.source)]
        ticks = percent_ticks(axis_limits, CORRELATION_TICK_STEPS[0])
        layout = CORRELATION_LAYOUTS[group]
        stacked = set_id in str(comparison.stacked_sets).split(MEASUREMENT_SEPARATOR)
        shape = (len(metadata), 1) if stacked else (1, len(metadata))
        figure, axes = plt.subplots(*shape, figsize=STACKED_SIZE if stacked else layout["size"], layout="constrained", sharex=True, sharey=True)
        for ax, reaction in zip(axes, metadata.itertuples()):
            density_ax = make_axes_locatable(ax).append_axes("right", size=f"{100 * DENSITY_WIDTH_RATIO}%", pad=DENSITY_GAP_IN, sharey=ax)
            points, dropped = comparison_points(set_id, comparison, reaction)
            x, y = points[comparison.x_position_id].to_numpy(), points[comparison.y_position_id].to_numpy()
            ax.plot(axis_limits, axis_limits, linestyle="--", color=REFERENCE_COLOR, linewidth=.4)
            ax.scatter(x, y, color=colors[reaction.tada_type], linewidths=0, **VARIANT_POINT_STYLES[comparison.variant_set])
            # Match R density's default nrd0 bandwidth, rather than scipy's Scott default.
            bandwidth = .9 * min(np.std(y, ddof=1), np.subtract(*np.percentile(y, [75, 25])) / 1.34) * len(y) ** -.2
            density_y = np.linspace(*axis_limits, DENSITY_POINTS)
            density = gaussian_kde(y, bw_method=bandwidth / np.std(y, ddof=1))(density_y)
            density_ax.fill_betweenx(density_y, 0, density / density.max(), color=colors[reaction.tada_type], alpha=DENSITY_ALPHA, linewidth=0)
            density_ax.axvline(0, color=REFERENCE_COLOR, linewidth=.4)
            density_ax.set_xlim(0, 1.05)
            density_ax.axis("off")
            # Only the correlation is drawn on the panel. How many variants it rests on belongs in
            # the caption, which can also say which ones are left out and why; the count per panel
            # is printed on every run and kept in the statistics table.
            ax.text(.04, .96, f"r = {np.corrcoef(x, y)[0, 1]:.2f}", transform=ax.transAxes, va="top")
            ax.set(title=titles[reaction.sample_id], xlim=axis_limits, ylim=axis_limits)
            ax.set_aspect("equal")
            # Ticks sit on powers of two, labelled as the percentage of wild-type editing, as every
            # other axis and colorbar in this analysis is.
            percent_labels = [f"{100 * 2. ** power:.3g}" for power in ticks]
            ax.set_xticks(ticks, labels=percent_labels)
            ax.set_yticks(ticks, labels=percent_labels)
            linear = np.corrcoef(100 * 2. ** x, 100 * 2. ** y)[0, 1]
            print(f"{set_id} {comparison.variant_set} {comparison.x_position_id} vs {comparison.y_position_id}, {titles[reaction.sample_id]}: r={np.corrcoef(x, y)[0, 1]:.10f} on log2 "
                  f"(r={linear:.10f} on percent), n={len(x)}, {dropped} variant(s) with no edited molecule left out, "
                  f"median={100 * 2. ** np.median(x):.10f}%, SD={np.std(x, ddof=1):.10f} log2")
            # Summarize each drawn axis on its own, as the percentage of wild-type editing the panel
            # labels and as the loss below wild type. The median and the geometric mean are taken on
            # the log₂ ratios the panel plots, which is where a fold change is a distance; the plain
            # mean of the percentages is kept beside them, since a few variants above wild type pull
            # it well above the other two and the difference is itself worth reading.
            for position_id, values in [(comparison.x_position_id, x), (comparison.y_position_id, y)]:
                percent = 100 * 2. ** values
                median_percent, geometric_mean_percent = 100 * 2. ** np.median(values), 100 * 2. ** values.mean()
                correlation_statistics.append(dict(
                    reaction_set=set_id, source=f"{comparison.source}_{set_id}", variant_set=comparison.variant_set,
                    comparison=f"{comparison.x_position_id} vs {comparison.y_position_id}", enzyme=titles[reaction.sample_id],
                    position_id=position_id, n_variants=len(values), n_zero_edit_left_out=dropped,
                    median_percent_of_wt=median_percent, median_decrease_percentage_points=100 - median_percent,
                    geometric_mean_percent_of_wt=geometric_mean_percent, geometric_mean_decrease_percentage_points=100 - geometric_mean_percent,
                    mean_percent_of_wt=percent.mean(), mean_decrease_percentage_points=100 - percent.mean(),
                    sd_log2=np.std(values, ddof=1), pearson_r_log2=np.corrcoef(x, y)[0, 1], pearson_r_percent=linear))
        # Side by side the panels share one y axis, drawn only at the left, and one x label, while each
        # keeps its own x ticks: the two enzymes are read against each other rather than one above the other.
        # Stacked they share the x axis the same way, drawn only under the lower panel.
        if stacked:
            for ax in axes[:-1]:
                ax.tick_params(bottom=False, labelbottom=False)
                ax.spines["bottom"].set_visible(False)
        else:
            for ax in axes[1:]:
                ax.tick_params(left=False, labelleft=False)
                ax.spines["left"].set_visible(False)
        # Only the panels that still draw x numbers can be measured for how far apart they stand.
        labelled = [axes[-1]] if stacked else list(axes)
        # A deep axis can bring its own numbers together faster than it brings the two panels together,
        # which the clearance check below does not see, so thin the lattice until one axis's own
        # neighbouring numbers stand apart. Only the designed panels reach that deep; the
        # single-substitution panels take the first step and are untouched.
        if comparison.variant_set == "designed":
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            for step in CORRELATION_TICK_STEPS:
                ticks = percent_ticks(axis_limits, step)
                percent_labels = [f"{100 * 2. ** power:.3g}" for power in ticks]
                for ax in axes:
                    ax.set_xticks(ticks, labels=percent_labels)
                    ax.set_yticks(ticks, labels=percent_labels)
                figure.canvas.draw()
                drawn = [[number.get_window_extent(renderer) for number in ax.get_xticklabels()] for ax in labelled]
                spacing = 72 * min(right.x0 - left.x1 for row in drawn for left, right in zip(row, row[1:])) / figure.dpi
                if spacing >= TICK_LABEL_GAP:
                    break
            require(spacing >= TICK_LABEL_GAP, f"One axis's own tick numbers stay within {TICK_LABEL_GAP} pt at every step in {comparison.source}_{set_id}")
            print(f"{set_id} {comparison.source}: labelled every {step} powers of two at "
                  + ", ".join(f"{100 * 2. ** power:g}%" for power in ticks)
                  + f", {spacing:.1f} pt between neighbouring numbers")
        # A stacked figure is a third as wide, so its x label breaks at its own comma as the y label does.
        shared_x = figure.supxlabel((stacked_labels if stacked else axis_labels)[comparison.x_position_id],
                                    multialignment="center")
        place_shared_ylabel(figure, axes[0], stacked_labels[comparison.y_position_id], gap=layout["ylabel_gap"])
        centre_shared_xlabel(figure, shared_x)
        figure.canvas.draw()
        label_box = shared_x.get_window_extent(figure.canvas.get_renderer())
        require(label_box.x0 >= 0 and label_box.x1 <= figure.bbox.x1, f"The x label does not fit the page in {comparison.source}_{set_id}")
        # Panels drawn this close must still keep their tick numbers apart, so measure the drawn
        # numbers and fail here rather than ship a figure whose axes read into each other.
        renderer = figure.canvas.get_renderer()
        if stacked:
            for upper, lower in zip(axes, axes[1:]):
                clearance = min(number.get_window_extent(renderer).y0 for number in upper.get_yticklabels()) - max(number.get_window_extent(renderer).y1 for number in lower.get_yticklabels())
                require(clearance >= 0, f"Stacked panels' tick numbers touch in {comparison.source}_{set_id}")
        else:
            for left, right in zip(axes, axes[1:]):
                clearance = min(number.get_window_extent(renderer).x0 for number in right.get_xticklabels()) - max(number.get_window_extent(renderer).x1 for number in left.get_xticklabels())
                require(clearance >= 0, f"Neighbouring panels' tick numbers touch in {comparison.source}_{set_id}")
        # The all-variant panels are read on their own rather than assembled beside the heat maps, so
        # they always carry a PNG companion for the README; the single-substitution sets keep theirs
        # with the layout they are assembled in.
        png = layout["png"] or comparison.variant_set == "designed"
        save_pdf(figure, f"{comparison.source}_{set_id}", transparent=VARIANT_TRANSPARENT[comparison.variant_set], png=png, compact=layout["compact"])

# One table for every correlation panel this family draws, so the numbers the figures carry can be
# read and checked without re-running the plotter. Written once, after every set, rather than
# appended per set, which would leave a partial table behind a failed run.
statistics = pd.DataFrame(correlation_statistics)
require(len(statistics) > 0 and statistics.notna().all().all() and not statistics.duplicated(["source", "enzyme", "position_id"]).any(), "Invalid correlation statistics")
require(np.isfinite(statistics.select_dtypes("number")).all().all(), "Non-finite correlation statistics")
os.makedirs(os.path.dirname(STATISTICS_PATH), exist_ok=True)
statistics.to_csv(STATISTICS_PATH, index=False)
print(f"Saved {STATISTICS_PATH}; {len(statistics)} panel axes over {statistics.source.nunique()} figures")
