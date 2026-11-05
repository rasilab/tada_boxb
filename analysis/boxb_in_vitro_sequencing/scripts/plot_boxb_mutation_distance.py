"""Distribution of boxB window effects by distance from wild type; sequences are not reaction replicates."""

import math
import numpy as np
from analysis_helpers import mutation_distance_data, reaction_set, reaction_set_ids, require, scale_groups
from plot_style import plt, palettes, save_pdf, place_shared_ylabel, centre_shared_xlabel, LABEL_GAP

# One layout per scale group; shorter editor panels compact the horizontal distribution figure.
LAYOUTS = {
    "lambdaN_tada": dict(columns=1, size=(1.55, 3.12), aspect=1, ylabel_gap=LABEL_GAP, compact=False, png=False),
    "editors": dict(columns=2, size=(2.20, 1.15), aspect=.8, ylabel_gap=1 / 72, compact=True, png=True),
    "time": dict(columns=1, size=(1.55, 3.12), aspect=1, ylabel_gap=LABEL_GAP, compact=False, png=True),
}
VIOLIN_WIDTH = .80
VIOLIN_ALPHA = .55
MEDIAN_WIDTH = .34
MEDIAN_LINEWIDTH = .9
CONTROL_POINT_SIZE = 1.6
CONTROL_JITTER = .20
CONTROL_SEED = 0
REFERENCE_COLOR = "#999999"
REFERENCE_WIDTH = .4
TICK_STEP = 2
LOWEST_TICK = -4
AXIS_PAD = .25
CONTROL_POSITION = 0
SOURCE = "boxb_mutation_distance"
TOLERANCE = 1e-9

labels, _, colors = palettes()
groups = scale_groups()
group_of = {set_id: group for group, members in groups.items() for set_id in members}
set_ids = reaction_set_ids()


def plotted(frame, column):
    """The values that can sit on a log axis.

    A sequence with no edited molecule has no logarithm, so it cannot be placed. It is left out of its
    violin and counted in the printed summary rather than imputed to the axis floor, which would invent
    a value for the very sequences the measurement could not reach.
    """
    values = frame[column].to_numpy()
    return values[np.isfinite(values)]


# Build every set before drawing any, so the sets of one scale group share a y axis and a violin can
# be compared across them. The shared axis is the widest any of them needs, so nothing is clipped.
# Effects are spaced in log₂ of the ratio to wild type, which gives equal fold changes up and down equal
# distance and puts wild type on the reference line, so one violin's spread can be compared with
# another's by eye; the ticks are labelled as the percentage of wild-type editing, on powers of two,
# exactly as the combinatorial colorbar and the correlation scatters are.
datasets = {set_id: mutation_distance_data(set_id) for set_id in set_ids}
require(len({measurement.position_id for _, _, measurement in datasets.values()}) == 1, "Every set must measure the same recorder site")
site = next(iter(datasets.values()))[2].axis_label
label = f"{site}, % of wild type"
# Beside the panels the label sits in a strip of its own, so it breaks at its comma as the scatters' does.
stacked_label = f"{site},\n% of wild type"
axis_limits, ticks = {}, {}
for group, members in groups.items():
    spread = [plotted(frame, column) for set_id in members for distances, controls, _ in [datasets[set_id]]
              for frame, column in [(distances, "log2_pooled_ratio_to_matched_wt"), (controls, "log2_ratio_to_pooled_wt")]]
    # The axis reaches the weakest sequence of the group rather than stopping at a floor: a violin cut
    # off against the axis reads as a complete distribution, and the lower tail is the measurement.
    # The labelled ticks stay where they are, on powers of two from 6.25% of wild type up, so these
    # panels still read against the rest of the set; the stretch below the lowest tick is unlabelled.
    combined = np.concatenate(spread)
    limits = (min(LOWEST_TICK, combined.min() - AXIS_PAD), combined.max() + AXIS_PAD)
    axis_limits[group] = limits
    ticks[group] = list(range(LOWEST_TICK, math.floor(limits[1]) + 1, TICK_STEP))
    print(f"{group}: the axis reaches the weakest of {len(combined)} effects, {100 * 2. ** combined.min():.2f}% of wild type, so no violin is cut off")
    require(len(ticks[group]) > 1, f"The shared axis for {group} carries no labelled fold change")
    print(f"{label}; {group} axis {100 * 2. ** limits[0]:.1f}–{100 * 2. ** limits[1]:.1f}% of wild type, "
          f"ticks at {', '.join(f'{100 * 2. ** power:g}' for power in ticks[group])}%")

# A mutant is matched to the wild type sharing its barcode, so its own ratio is one by construction and
# carries no spread. Position zero therefore holds the wild-type controls against their pooled mean,
# which is the baseline every mutant ratio is taken against, rather than a degenerate point at zero.
distance_values = sorted({int(value) for frame, _, _ in datasets.values() for value in frame.n_mutations if value > 0})
positions = [CONTROL_POSITION] + distance_values
for set_id in set_ids:
    group = group_of[set_id]
    frame, controls, _ = datasets[set_id]
    # One panel per reaction rather than per enzyme, since a set may draw one enzyme twice.
    enzymes = reaction_set(set_id)
    titles = {row.sample_id: row.panel_title if isinstance(row.panel_title, str) else labels[row.tada_type] for row in enzymes.itertuples()}
    layout = LAYOUTS[group]
    columns = layout["columns"]
    figure, axes = plt.subplots(len(enzymes) // columns, columns, figsize=layout["size"], layout="constrained", sharex=True, sharey=True)
    for ax, enzyme in zip(axes, enzymes.itertuples()):
        rows = frame.loc[frame.sample_id.eq(enzyme.sample_id)]
        control_values = plotted(controls.loc[controls.sample_id.eq(enzyme.sample_id)], "log2_ratio_to_pooled_wt")
        by_distance = [rows.loc[rows.n_mutations.eq(distance)] for distance in distance_values]
        violins = [control_values] + [plotted(part, "log2_pooled_ratio_to_matched_wt") for part in by_distance]
        unplaceable = sum(len(part) for part in by_distance) - sum(len(values) for values in violins[1:])
        require(all(len(values) > 1 for values in violins), f"A position carries too few sequences to shape a violin in {set_id}")
        ax.axhline(0, linestyle="--", color=REFERENCE_COLOR, linewidth=REFERENCE_WIDTH)
        parts = ax.violinplot(violins, positions=positions, widths=VIOLIN_WIDTH, showextrema=False, showmedians=False)
        for body in parts["bodies"]:
            body.set(facecolor=colors[enzyme.tada_type], edgecolor="none", alpha=VIOLIN_ALPHA)
        # The controls are few enough to show outright, so the reader is not asked to trust a density
        # smoothed over 22 points the way the hundreds behind each mutant violin can be.
        offsets = np.random.default_rng(CONTROL_SEED).uniform(-CONTROL_JITTER, CONTROL_JITTER, len(control_values))
        ax.plot(CONTROL_POSITION + offsets, control_values, linestyle="none", marker="o", markersize=CONTROL_POINT_SIZE,
                markerfacecolor=colors[enzyme.tada_type], markeredgecolor="none")
        # The median is the one summary the violin does not already show, so mark it and nothing else.
        for position, values in zip(positions, violins):
            median = np.median(values)
            ax.plot([position - MEDIAN_WIDTH, position + MEDIAN_WIDTH], [median, median], color=colors[enzyme.tada_type], linewidth=MEDIAN_LINEWIDTH, solid_capstyle="butt")
        ax.set(title=titles[enzyme.sample_id], ylim=axis_limits[group], yticks=ticks[group], yticklabels=[f"{100 * 2. ** power:g}" for power in ticks[group]],
               xlim=(min(positions) - .5 - AXIS_PAD, max(positions) + .5 + AXIS_PAD), xticks=positions)
        ax.set_box_aspect(layout["aspect"])
        print(f"{set_id} {titles[enzyme.sample_id]}: " + ", ".join(
            f"{'WT controls' if position == CONTROL_POSITION else f'{position} subs'} n={len(values)} "
            f"median={100 * 2. ** np.median(values):.1f}% IQR {100 * 2. ** np.percentile(values, 25):.1f}–{100 * 2. ** np.percentile(values, 75):.1f}%"
            for position, values in zip(positions, violins))
            + (f"; {unplaceable} sequence(s) with no edited molecule left out of their violin" if unplaceable else ""))
    if columns == 1:
        # The bottom panel carries the shared x axis; the top panel drops its own.
        axes[0].tick_params(bottom=False, labelbottom=False)
        axes[0].spines["bottom"].set_visible(False)
        # The label breaks onto two lines, since on a panel this narrow one line runs off the page.
        axes[-1].set_xlabel("Substitutions from\nwild type")
        drawn = axes[-1].xaxis.label
        place_shared_ylabel(figure, axes[0], label, gap=layout["ylabel_gap"])
    else:
        # Side by side the panels share one y axis, drawn only at the left, and one x label under both,
        # while each keeps its own x ticks: the two editors are read against each other.
        axes[-1].tick_params(left=False, labelleft=False)
        axes[-1].spines["left"].set_visible(False)
        drawn = figure.supxlabel("Substitutions from wild type")
        place_shared_ylabel(figure, axes[0], stacked_label, gap=layout["ylabel_gap"])
        centre_shared_xlabel(figure, drawn)
    box = drawn.get_window_extent(figure.canvas.get_renderer())
    require(box.x0 >= 0 and box.x1 <= figure.bbox.x1, f"The x label runs off the page in {set_id}")
    save_pdf(figure, f"{SOURCE}_{set_id}", transparent=True, png=layout["png"], compact=layout["compact"])
