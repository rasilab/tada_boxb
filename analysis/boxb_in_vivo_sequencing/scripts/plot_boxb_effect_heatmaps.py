"""Aligned single-mutant heat map and downward mean-editing-loss logo.

Edited and callable counts are pooled over matching constructs separately for each site.
The resulting editing percentages are averaged with equal site weights, retaining observed zeros.
WT uses its measured editing percentage on the same linear scale as every single mutant.
Each position's four-base mean is expressed as a percentage reduction from the highest position mean.
Downward stacks divide that loss among upright letters in proportion to their editing percentages.
"""

import os

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import ScaledTranslation

from analysis_helpers import (DNA_BASES, annotation, boxb_wildtype, comparison_samples,
                              designed_boxb_sequences, excluded_constructs,
                              require, sequence_effects, substitution_measurements)
from plot_style import logo_letter, palettes, plt, save_figure

CELL_HEIGHT = .14
FIGURE_WIDTH = 5.1  # Base layout; extend only the left margin for the single-line annotation labels.
FIGURE_HEIGHT = 2.0
ANNOTATION_HEIGHT = .62
LOGO_HEIGHT = .60
COLOR_TICK_BINS = 3
TICK_INTERVALS = [1, 2, 2.5, 5, 10]
TICK_PAD = 3
YLABEL_PAD = 1
LOGO_LABEL_SHIFT_PT = 5
HIGH_COLOR = "white"
DARK_TEXT = "#333333"
LIGHT_TEXT = "white"
LUMINANCE_WEIGHTS = [.2126, .7152, .0722]
TEXT_CONTRAST_THRESHOLD = .5
TILE_FONT = 7
LOGO_WIDTH = .85
LOGO_LIMITS = (-100, 0)
LOGO_TICKS = [-50, 0]
LOGO_YLABEL = "Editing loss\n(% of max)"
ZERO_LINE_WIDTH = .4
ZERO_LINE_COLOR = "#777777"
HEAT_LABEL_X = 1.02
TABLE_DIRECTORY = "analysis/boxb_in_vivo_sequencing/data/figure_data"
LAYOUT_PADDING = .02
LAYOUT_HSPACE = .02
ALIGNMENT_TOLERANCE = 1e-8
# Preserve the in vitro annotation rows and span geometry with tighter vertical spacing.
BASE_OFFSET_PT = 8
POSITION_OFFSET_PT = 24
POSITION_LABEL_SHIFT_PT = 3
STRUCTURE_OFFSET_PT = 38
ROW_LABEL_X = -.9
STRUCTURE_BAR_PT = 4
STRUCTURE_GAP = .25
STRUCTURE_EDGE_INSET = .1
STRUCTURE_LINE_WIDTH = .4
ANNOTATION_COLOR = "black"

labels, _, colors = palettes()
selections = annotation("boxb_effect_heatmaps")
require(selections.notna().all().all() and not selections.source.duplicated().any(), "Invalid heat-map selections")
require(excluded_constructs().empty, "Flagged placements require an explicit with/without comparison")
regions, wildtype_rna, wildtype = boxb_wildtype()
designed = designed_boxb_sequences()
singles = designed.loc[designed.n_mutations.le(1)]
grid = pd.DataFrame([dict(boxb_position=index + 1, nucleotide=base.replace("T", "U"),
                          is_wildtype=base == reference, hairpin=wildtype[:index] + base + wildtype[index + 1:])
                     for index, reference in enumerate(wildtype) for base in DNA_BASES])
require(set(grid.hairpin) == set(singles.hairpin), "The substitution grid and designed library disagree")
rna_bases = [base.replace("T", "U") for base in DNA_BASES]
prepared = {}

for (set_id, sample_id), definitions in selections.groupby(["comparison_set", "sample_id"], sort=False):
    selected = comparison_samples(set_id)
    selected = selected.loc[selected.sample_id.eq(sample_id)]
    require(len(selected) == 1, "A heat map must select one annotated library")
    sites = list(dict.fromkeys(";".join(definitions.sites).split(";")))
    measurements = substitution_measurements(set_id).set_index("position_id")
    require(set(sites) <= set(measurements.index), "An unannotated site was requested")
    measurements = measurements.loc[sites].reset_index()
    require(measurements.sites.eq(measurements.position_id).all(), "Average marginal site effects, not joint events")
    effects, constructs = sequence_effects(selected, measurements)
    effects = effects.merge(singles[["hairpin"]], on="hairpin", validate="many_to_one")
    require(len(effects) == len(selected) * len(sites) * len(singles)
            and not effects.duplicated(["sample_id", "position_id", "hairpin"]).any(),
            "A library is missing a designed single-mutant or WT measurement")
    require(effects.pooled_umis.gt(0).all() and effects.pooled_fraction_edited.between(0, 1).all(),
            "A site lacks callable coverage or a valid editing fraction")

    for definition in definitions.itertuples():
        members = definition.sites.split(";")
        require(len(members) == len(set(members)), "An average repeats a site")
        chosen = effects.loc[effects.position_id.isin(members)]
        average = chosen.groupby(["sample_id", "hairpin"], as_index=False).agg(
            mean_fraction_edited=("pooled_fraction_edited", "mean"), n_sites=("position_id", "size"))
        require(average.n_sites.eq(len(members)).all(), "Every effect must include every requested site")
        tiles = average.merge(grid, on="hairpin", validate="many_to_many")
        require(len(tiles) == len(selected) * len(grid), "Incomplete single-mutant heat map")
        tiles["mean_editing_percent"] = 100 * tiles.mean_fraction_edited
        prepared[definition.source] = tiles

    heat_height = len(rna_bases) * CELL_HEIGHT
    for definition in definitions.itertuples():
        tiles = prepared[definition.source]
        position_totals = tiles.groupby(["sample_id", "boxb_position"]).mean_editing_percent.transform("sum")
        require(position_totals.gt(0).all(), "A zero-edit position cannot have a normalized logo")
        tiles["relative_editing"] = tiles.mean_editing_percent / position_totals
        tiles["position_mean_editing_percent"] = tiles.groupby(
            ["sample_id", "boxb_position"]).mean_editing_percent.transform("mean")
        tiles["maximum_position_mean_editing_percent"] = tiles.groupby(
            "sample_id").position_mean_editing_percent.transform("max")
        tiles["position_editing_change_percent"] = 100 * (
            tiles.position_mean_editing_percent / tiles.maximum_position_mean_editing_percent - 1)
        tiles["logo_letter_change_percent"] = tiles.relative_editing * tiles.position_editing_change_percent
        require(np.allclose(tiles.groupby(["sample_id", "boxb_position"]).relative_editing.sum(), 1),
                "Logo columns do not sum to one")
        position_changes = tiles.groupby(["sample_id", "boxb_position"])
        require(tiles.position_editing_change_percent.between(-100, 0).all()
                and np.allclose(tiles.groupby("sample_id").position_editing_change_percent.max(), 0)
                and np.allclose(position_changes.logo_letter_change_percent.sum(),
                                position_changes.position_editing_change_percent.first()),
                "Logo stacks must sum to the percentage change from each library's maximum position mean")
        os.makedirs(TABLE_DIRECTORY, exist_ok=True)
        tiles.to_csv(f"{TABLE_DIRECTORY}/{definition.source}.csv", index=False)
        figure = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT), layout="constrained")
        figure.get_layout_engine().set(w_pad=LAYOUT_PADDING, h_pad=LAYOUT_PADDING, hspace=LAYOUT_HSPACE)
        layout = figure.add_gridspec(3, 1, height_ratios=[heat_height, LOGO_HEIGHT, ANNOTATION_HEIGHT])
        for row, library in enumerate(selected.itertuples()):
            ax = figure.add_subplot(layout[0, 0])
            logo_ax = figure.add_subplot(layout[1, 0], sharex=ax)
            annotation_ax = figure.add_subplot(layout[2, 0], sharex=ax)
            annotation_ax.set_ylim(0, 1)
            annotation_ax.set_axis_off()
            data = tiles.loc[tiles.sample_id.eq(library.sample_id)]
            ticks = MaxNLocator(nbins=COLOR_TICK_BINS, steps=TICK_INTERVALS).tick_values(0, data.mean_editing_percent.max())
            limits = (0, ticks[-1])
            require(limits[1] > limits[0], "The effects have no color range")
            norm = Normalize(*limits)
            percent = data.pivot(index="nucleotide", columns="boxb_position", values="mean_editing_percent").loc[rna_bases].to_numpy()
            wt = data.pivot(index="nucleotide", columns="boxb_position", values="is_wildtype").loc[rna_bases].to_numpy()
            zero = percent == 0
            cmap = LinearSegmentedColormap.from_list("enzyme", [colors[library.enzyme], HIGH_COLOR])
            ax.pcolormesh(np.arange(len(wildtype) + 1) - .5, np.arange(len(rna_bases) + 1) - .5,
                          percent, cmap=cmap, norm=norm, rasterized=False)
            for (y, x), value in np.ndenumerate(percent):
                background = cmap(norm(value))[:3]
                text_color = LIGHT_TEXT if np.dot(background, LUMINANCE_WEIGHTS) < TEXT_CONTRAST_THRESHOLD else DARK_TEXT
                ax.text(x, y, f"{value:.1f}", ha="center", va="center", color=text_color, fontsize=TILE_FONT)
            ax.set(xlim=(-.5, len(wildtype) - .5), ylim=(len(rna_bases) - .5, -.5),
                   yticks=range(len(rna_bases)), yticklabels=rna_bases, ylabel="Mutant\nbase")
            ax.tick_params(length=0, pad=TICK_PAD)
            ax.tick_params(axis="x", labelbottom=False)
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.text(HEAT_LABEL_X, .5, definition.label.replace(" editing", "\nediting"),
                    rotation=270, va="center", ha="left", multialignment="center", transform=ax.transAxes)
            relative = data.pivot(index="nucleotide", columns="boxb_position", values="relative_editing").loc[rna_bases]
            changes = data.groupby("boxb_position").position_editing_change_percent.first()
            reference = data.maximum_position_mean_editing_percent.iloc[0]
            for position in relative.columns:
                bottom = 0
                # Largest letters sit nearest zero; positive glyph heights keep every letter upright.
                for base, fraction in relative[position].sort_values(ascending=False, kind="stable").items():
                    height = -changes.loc[position] * fraction
                    bottom -= height
                    if height > 0:
                        logo_letter(logo_ax, base, position - 1, bottom, LOGO_WIDTH, height)
                require(np.isclose(bottom, changes.loc[position]), "A drawn logo does not match its editing reduction")
            logo_ax.axhline(0, color=ZERO_LINE_COLOR, linewidth=ZERO_LINE_WIDTH, clip_on=False)
            # Trim only the visible axis to the deepest stack, preserving the logo's scale and space.
            logo_ax.set(ylim=LOGO_LIMITS, yticks=LOGO_TICKS, ylabel=LOGO_YLABEL)
            logo_ax.spines["left"].set_bounds(changes.min(), 0)
            logo_ax.tick_params(axis="x", length=0, labelbottom=False)
            logo_ax.spines["bottom"].set_visible(False)
            row_labels = []
            for offset, label_shift, row_label, row_values in [
                    (BASE_OFFSET_PT, 0, "Wild-type\nbase", wildtype_rna),
                    (POSITION_OFFSET_PT, POSITION_LABEL_SHIFT_PT, "boxB\nposition", range(1, len(wildtype) + 1))]:
                transform = annotation_ax.get_xaxis_transform() + ScaledTranslation(
                    0, -offset / 72, figure.dpi_scale_trans)
                for x, value in enumerate(row_values):
                    annotation_ax.text(x, 1, str(value), ha="center", color=ANNOTATION_COLOR, transform=transform)
                label_transform = transform + ScaledTranslation(0, -label_shift / 72, figure.dpi_scale_trans)
                text = annotation_ax.text(ROW_LABEL_X, 1, row_label, ha="right", va="center", multialignment="right",
                                          color=ANNOTATION_COLOR, transform=label_transform)
                row_labels.append((text, transform))
            structure_transform = annotation_ax.get_xaxis_transform() + ScaledTranslation(
                0, -STRUCTURE_OFFSET_PT / 72, figure.dpi_scale_trans)
            structure_labels = []
            for region in regions.itertuples():
                left = region.start - 1.5 + STRUCTURE_EDGE_INSET
                right = region.end - .5 - STRUCTURE_EDGE_INSET
                text = annotation_ax.text((left + right) / 2, 1, region.label, ha="center", va="center",
                                          color=ANNOTATION_COLOR, transform=structure_transform)
                structure_labels.append((left, right, text))
            print(f"{definition.source}, {labels[library.enzyme]}: {len(singles) - 1} mutants; "
                  f"{int(zero.sum())} observed zeros; "
                  f"range {percent[~wt].min():.3f}-{percent[~wt].max():.3f}% editing; "
                  f"WT {percent[wt][0]:.3f}%; linear color limits {limits}")
            print(f"Logo zero reference: position(s) {changes.index[changes.eq(0)].tolist()}, "
                  f"mean editing {reference:.6f}%; largest reduction {changes.min():.3f}% "
                  f"at position {changes.idxmin()}")
        labeled_axes = [ax, logo_ax]
        for panel_ax in labeled_axes:
            panel_ax.yaxis.labelpad = YLABEL_PAD
            panel_ax.yaxis.label.set_verticalalignment("bottom")
            panel_ax.yaxis.label.set_multialignment("center")
        figure.align_ylabels(labeled_axes)
        figure.canvas.draw()
        figure.set_layout_engine("none")
        renderer = figure.canvas.get_renderer()
        # Keep the existing panel geometry and horizontal y-label offsets, then extend the page left.
        logo_ax.yaxis.label.set_transform(logo_ax.yaxis.label.get_transform() + ScaledTranslation(
            0, LOGO_LABEL_SHIFT_PT / 72, figure.dpi_scale_trans))
        for text, transform in row_labels:
            text.set_text(text.get_text().replace("\n", " "))
            text.set_verticalalignment("baseline")
            text.set_transform(transform)
        left_edge = min(text.get_window_extent(renderer).x0 for text, _ in row_labels) / figure.dpi
        extra_left = max(0, LAYOUT_PADDING - left_edge)
        width, height = figure.get_size_inches()
        panel_positions = [panel.get_position().bounds for panel in figure.axes]
        figure.set_size_inches(width + extra_left, height)
        for panel, (x, y, w, h) in zip(figure.axes, panel_positions):
            panel.set_position(((x * width + extra_left) / (width + extra_left), y,
                                w * width / (width + extra_left), h))
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        print(f"Figure canvas {width + extra_left:.3f} × {height:.3f} inches; "
              f"added {extra_left:.3f} inches at left; raised logo label {LOGO_LABEL_SHIFT_PT} points")
        label_boxes = [panel_ax.yaxis.label.get_window_extent(renderer) for panel_ax in labeled_axes]
        require(np.allclose([box.x1 for box in label_boxes], label_boxes[0].x1),
                "The rotated y-label bottom edges are not aligned")
        require(all(upper.y0 > lower.y1 for upper, lower in zip(label_boxes, label_boxes[1:])),
                "Adjacent y-axis labels overlap")
        for left, right, text in structure_labels:
            box = text.get_window_extent(renderer)
            start, end = (annotation_ax.transData.inverted().transform((x, 0))[0] for x in [box.x0, box.x1])
            require(left < start - STRUCTURE_GAP and end + STRUCTURE_GAP < right,
                    f"Region label does not fit its span: {text.get_text()}")
            for span in [[left, start - STRUCTURE_GAP], [end + STRUCTURE_GAP, right]]:
                annotation_ax.plot(span, [1, 1], color=ANNOTATION_COLOR, linewidth=STRUCTURE_LINE_WIDTH,
                                   transform=structure_transform, clip_on=False)
            annotation_ax.plot([left, right], [1, 1], linestyle="none", marker="|", markersize=STRUCTURE_BAR_PT,
                               markeredgewidth=STRUCTURE_LINE_WIDTH, color=ANNOTATION_COLOR,
                               transform=structure_transform, clip_on=False)
        require(all(np.allclose([other.get_position().x0, other.get_position().x1],
                                [ax.get_position().x0, ax.get_position().x1], atol=ALIGNMENT_TOLERANCE)
                    for other in [annotation_ax, logo_ax]), "The position axes are not aligned")
        for ax in figure.axes:
            texts = list(ax.texts)
            if ax.axison:
                texts += [ax.xaxis.label, ax.yaxis.label, ax.title, *ax.get_xticklabels(), *ax.get_yticklabels()]
            for text in texts:
                if text.get_visible() and text.get_text():
                    box = text.get_window_extent(renderer)
                    require(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width
                            and box.y1 <= figure.bbox.height,
                            f"Figure label {text.get_text()!r} extends outside the page: {box.bounds}")
        save_figure(figure, definition.source)
    print("No in vivo placements are excluded; with/without-exclusion estimates are identical.")
