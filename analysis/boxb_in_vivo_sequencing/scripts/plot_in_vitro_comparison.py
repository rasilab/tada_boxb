"""Per-variant agreement between purified-enzyme and cellular editing of the boxB library.

Both axes carry the same quantity, built the same way on each side: a designed boxB sequence's
editing as a ratio to its matched wild type, taken as the plain mean of the A8 and A2 log2 ratios,
which is the site average the in vitro analysis defines. Normalising to matched wild type removes
the activity difference between a two-hour reaction and a cellular library and leaves the sequence
determinants, which is what the two experiments are being compared on.

The axes are the log2 ratio, labelled as the percentage of wild-type editing a reader can act on,
as every other variant correlation in this manuscript is drawn. A sequence with no edited molecule
has no logarithm on the axis it would be drawn against; it is left out and counted here rather than
carried in at a pseudocount, and the count is reported with the figure.

The pooled cloud is drawn first, as the agreement over the whole designed library, and then split by
the two kinds of randomized window the library holds. A window overlapping the loop holds boxB's own
adenosines, so an edit inside it changes the genotype its molecule is counted as; a stem window holds
none, and an edit there is recorded as recorder self-editing instead. The split panels show what the
pooled coefficient is composed of.
"""

import math
import os

import numpy as np
from scipy.stats import pearsonr, spearmanr

from analysis_helpers import (annotation, boxb_wildtype, comparison_samples, designed_boxb_sequences, designed_grid,
                              excluded_constructs, require, sequence_effects,
                              substitution_measurements)
from in_vitro_reuse import helpers as in_vitro
from plot_style import centre_shared_xlabel, palettes, place_shared_ylabel, plt, save_figure

# Square panels side by side, 0.925 inch each, since this is read as a supplementary panel rather than
# beside the main figures. The canvas is sized from that panel: the label and tick margins it leaves
# are fixed in inches, so the width is those margins plus one panel per set drawn, and the height
# those margins plus one.
FIGURE_SIZE = (4.49, 1.56)
# Every panel draws the same marker, so a reader compares the panels on where their points sit rather
# than on how they are drawn. The size and opacity are a compromise: large and opaque enough that the
# fifty-odd points of the sparse single-substitution panel read clearly, small and light enough that
# the hundreds of points of a dense panel stay readable where they pile up.
POINT_STYLE = dict(s=3.5, alpha=.7)
# Each panel names the sequences it draws, in the order the annotation lists them: a title, the number
# of substitutions a sequence must carry, and whether its randomized window overlaps the loop. A
# selection left unset is not applied, so `designed` draws every designed sequence. Wild type carries
# no randomized window of its own and so falls out of the window panels.
VARIANT_SETS = {
    "designed": dict(title="All variants", n_mutations=None, loop_window=None),
    "single": dict(title="Single substitutions", n_mutations=1, loop_window=None),
    "stem": dict(title="Stem windows", n_mutations=None, loop_window=False),
    "loop": dict(title="Loop windows", n_mutations=None, loop_window=True),
}
SET_SEPARATOR = ";"
# The region whose adenosines a randomized window either does or does not contain.
LOOP_REGION = "loop"
# Half of wild-type editing, on the log2 axis. It is the level the summaries report a variant as
# falling below, and the level under which the reaction keeps almost no range of its own, so one
# constant serves both: a median alone cannot say how many variants lost editing, and the reaction's
# floor is what decides how much of a panel's spread is there to be correlated at all.
HALF_WILD_TYPE_LOG2 = -1.
# Ticks sit on powers of two and are labelled as percentages of wild type, with 100% always labelled.
TICK_STEP = 2
LOWEST_TICK = -4
AXIS_PADDING = .3
REFERENCE_COLOR = "#999999"
REFERENCE_WIDTH = .4
CORRELATION_POSITION = (.04, .96)
# The count sits in the opposite corner from the coefficient: the dense panel's cloud runs up the
# diagonal, so a two-line label in one corner lands on the points at its lower end.
COUNT_POSITION = (.96, .04)
SITE_SEPARATOR = ";"
SOURCE_DATA_DIRECTORY = "source_data"

_, _, colors = palettes()
selections = annotation("boxb_in_vitro_comparison")
require(selections.notna().all().all() and not selections.source.duplicated().any(),
        "Invalid or repeated in vitro comparison selections")
require(excluded_constructs().empty, "Flagged placements require an explicit with/without comparison")
designed = designed_boxb_sequences()

# Which of the two window kinds each designed mutant comes from, taken from the design and the
# annotated structure rather than from any count. A window spanning any loop position is a loop
# window, and the window a mutant belongs to is the one its substitutions were drawn from, not the
# region those substitutions happen to land in: it is the window that fixes which positions were
# randomized, and so whether boxB's own adenosines sit inside the counted insert.
regions, _, wildtype = boxb_wildtype()
loop_region = regions.loc[regions.region.eq(LOOP_REGION)]
require(len(loop_region) == 1, "The structure annotation must name one loop")
cells = designed_grid()
cells["loop_window"] = (cells.window_start.le(loop_region.iloc[0].end)
                        & cells.window_end.ge(loop_region.iloc[0].start))
window_kinds = cells.loc[cells.hairpin.ne(wildtype)].groupby("hairpin").loop_window.agg(["min", "max"])
require(window_kinds["min"].eq(window_kinds["max"]).all(),
        "A designed mutant comes from both a stem and a loop window")
loop_window = window_kinds["min"]


def site_mean_log2(effects, keys):
    """The in vitro site average applied unchanged: the plain mean of the log2 ratios.

    The sites are read from the same molecules, so their counts are never added, and a sequence
    unedited at a site carries minus infinity into the mean rather than a finite stand-in.
    """
    averaged = effects.groupby(keys, as_index=False).agg(
        log2_ratio=("log2_pooled_ratio_to_matched_wt", "mean"), n_sites=("position_id", "size"))
    return averaged


def percent_ticks(limits, step):
    """Powers-of-two ticks anchored on wild type, so 100% always carries a label."""
    lower, upper = limits
    return list(range(step * math.ceil(min(LOWEST_TICK, lower) / step), step * math.floor(upper / step) + 1, step))


for selection in selections.itertuples():
    selected = comparison_samples(selection.comparison_set)
    selected = selected.loc[selected.sample_id.eq(selection.sample_id)]
    require(len(selected) == 1, "A comparison must select one cellular library")
    sites = selection.sites.split(SITE_SEPARATOR)
    measurements = substitution_measurements(selection.comparison_set).set_index("position_id")
    require(len(sites) == len(set(sites)) and set(sites) <= set(measurements.index), "Invalid comparison sites")
    measurements = measurements.loc[sites].reset_index()
    require(measurements.sites.eq(measurements.position_id).all(), "Average marginal sites, not joint events")

    effects, constructs = sequence_effects(selected, measurements)
    effects = effects.merge(designed[["hairpin", "n_mutations"]], on=["hairpin", "n_mutations"], validate="many_to_one")
    require(len(effects) == len(designed) * len(sites)
            and not effects.duplicated(["hairpin", "position_id"]).any(), "Missing designed sequence measurement")
    require(effects.pooled_umis.gt(0).all() and effects.pooled_matched_wt_edited_umis.gt(0).all(),
            "A sequence lacks callable coverage or a matched wild type with edited molecules")
    in_cells = site_mean_log2(effects, ["hairpin", "n_mutations"])
    require(in_cells.n_sites.eq(len(sites)).all(), "Every sequence must include every selected site")

    measurement = in_vitro.selected_measurement(selection.in_vitro_measurement)
    reactions = in_vitro.reaction_set(selection.in_vitro_reaction_set)
    vitro_effects, _ = in_vitro.sequence_effects(reactions, measurement)
    vitro = vitro_effects.loc[vitro_effects.sample_id.eq(selection.in_vitro_sample_id)
                              & vitro_effects.position_id.eq(measurement.iloc[0].position_id)]
    require(len(vitro) == len(designed) and not vitro.hairpin.duplicated().any(),
            "The purified-enzyme reaction is missing a designed sequence")

    points = in_cells.merge(vitro[["hairpin", "log2_pooled_ratio_to_matched_wt"]].rename(
        columns={"log2_pooled_ratio_to_matched_wt": "in_vitro_log2"}), on="hairpin", validate="one_to_one")
    require(len(points) == len(designed), "The two experiments do not share every designed sequence")
    points["loop_window"] = points.hairpin.map(loop_window)
    require(points.loop_window.isna().sum() == 1 and points.loc[points.loop_window.isna()].n_mutations.eq(0).all(),
            "Only wild type carries no randomized window of its own")
    drawn = points.loc[np.isfinite(points.in_vitro_log2) & np.isfinite(points.log2_ratio)]
    exported = drawn.rename(columns={"hairpin": "boxb_sequence_dna", "log2_ratio": "in_cells_log2_ratio_to_wildtype",
                                     "in_vitro_log2": "in_vitro_log2_ratio_to_wildtype"}).copy()
    exported["loop_window"] = exported.loop_window.fillna("wild type")
    exported["in_cells_percent_of_wildtype"] = 100 * 2 ** exported.in_cells_log2_ratio_to_wildtype
    exported["in_vitro_percent_of_wildtype"] = 100 * 2 ** exported.in_vitro_log2_ratio_to_wildtype
    os.makedirs(SOURCE_DATA_DIRECTORY, exist_ok=True)
    exported.to_csv(f"{SOURCE_DATA_DIRECTORY}/{selection.source_data_panel}.csv", index=False)
    dropped = len(points) - len(drawn)
    x, y = drawn.in_vitro_log2.to_numpy(), drawn.log2_ratio.to_numpy()
    require(np.std(x) > 0 and np.std(y) > 0, "A constant measurement has no correlation")

    correlation = pearsonr(x, y)
    print(f"{selection.source}: n={len(drawn)}/{len(designed)} designed sequences, "
          f"{dropped} with no edited molecule left out; "
          f"Pearson r={correlation.statistic:.6f} on log2, p={correlation.pvalue:.3g}; "
          f"Spearman rho={spearmanr(x, y).statistic:.6f}")
    for label, mask in [("single substitutions", drawn.n_mutations.eq(1)),
                        ("two or more substitutions", drawn.n_mutations.ge(2)),
                        ("stem-window mutants", drawn.loop_window.eq(False)),
                        ("loop-window mutants", drawn.loop_window.eq(True))]:
        subset = drawn.loc[mask]
        if len(subset) > 2:
            print(f"  {label}: n={len(subset)}, "
                  f"Pearson r={pearsonr(subset.in_vitro_log2, subset.log2_ratio).statistic:.6f} on log2; "
                  f"median in vitro={100 * 2. ** subset.in_vitro_log2.median():.3f}% of WT, "
                  f"median in cells={100 * 2. ** subset.log2_ratio.median():.3f}% of WT")
            # The median says where the middle of a set sits, not how many of its variants lost
            # editing, and the two answer differently wherever a set mixes populations, as the four
            # stem windows do. Both are reported so neither can be read as the other.
            print(f"    below half of wild type: in vitro {subset.in_vitro_log2.lt(HALF_WILD_TYPE_LOG2).mean():.3f}, "
                  f"in cells {subset.log2_ratio.lt(HALF_WILD_TYPE_LOG2).mean():.3f}; "
                  f"SD in vitro={subset.in_vitro_log2.std(ddof=1):.3f}, in cells={subset.log2_ratio.std(ddof=1):.3f}; "
                  f"mean displacement={subset.log2_ratio.sub(subset.in_vitro_log2).mean():.3f} log2")

    # How much range each axis still carries where the reaction has stopped responding. The reaction
    # never falls near zero, so a set whose variants sit below half of wild type offers almost no in
    # vitro spread to correlate against however well the cells resolve them, and a coefficient read
    # over such a set measures the few variants that still edit rather than the whole panel.
    for label, mask in [("stem-window mutants", drawn.loop_window.eq(False)),
                        ("loop-window mutants", drawn.loop_window.eq(True))]:
        subset = drawn.loc[mask]
        print(f"  {label}: lowest in vitro {100 * 2. ** subset.in_vitro_log2.min():.3f}% of WT, "
              f"lowest in cells {100 * 2. ** subset.log2_ratio.min():.3f}% of WT")
        for part_label, keep in [("above half of wild type in vitro", subset.in_vitro_log2.ge(HALF_WILD_TYPE_LOG2)),
                                 ("below half of wild type in vitro", subset.in_vitro_log2.lt(HALF_WILD_TYPE_LOG2))]:
            part = subset.loc[keep]
            if len(part) > 2:
                print(f"    {part_label}: n={len(part)}, "
                      f"Pearson r={pearsonr(part.in_vitro_log2, part.log2_ratio).statistic:.6f} on log2; "
                      f"SD in vitro={part.in_vitro_log2.std(ddof=1):.3f}, "
                      f"in cells={part.log2_ratio.std(ddof=1):.3f}")

    # Every panel takes one pair of axes, spanning every drawn point, so the correlations are compared
    # on one scale rather than each rescaled to its own spread.
    library = selected.iloc[0]
    values = np.concatenate([x, y])
    limits = (values.min() - AXIS_PADDING, values.max() + AXIS_PADDING)
    ticks = percent_ticks(limits, TICK_STEP)
    tick_labels = [f"{100 * 2. ** power:.3g}" for power in ticks]
    require(0 in ticks, "Wild type must carry a tick label")
    variant_sets = selection.variant_sets.split(SET_SEPARATOR)
    require(len(variant_sets) == len(set(variant_sets)) and set(variant_sets) <= set(VARIANT_SETS),
            "A comparison must name distinct known variant sets")
    figure, axes = plt.subplots(1, len(variant_sets), figsize=FIGURE_SIZE, layout="constrained",
                                sharex=True, sharey=True, squeeze=False)
    for ax, variant_set in zip(axes.flat, variant_sets):
        panel_set = VARIANT_SETS[variant_set]
        title = panel_set["title"]
        selected_points = points
        if panel_set["n_mutations"] is not None:
            selected_points = selected_points.loc[selected_points.n_mutations.eq(panel_set["n_mutations"])]
        if panel_set["loop_window"] is not None:
            selected_points = selected_points.loc[selected_points.loop_window.eq(panel_set["loop_window"])]
        panel = drawn.loc[drawn.index.isin(selected_points.index)]
        require(len(panel) > 2, f"Too few sequences to correlate in the {variant_set} panel")
        panel_x, panel_y = panel.in_vitro_log2.to_numpy(), panel.log2_ratio.to_numpy()
        panel_correlation = pearsonr(panel_x, panel_y)
        left_out = len(selected_points) - len(panel)
        print(f"  {title}: n={len(panel)}/{len(selected_points)}, "
              f"Pearson r={panel_correlation.statistic:.6f} on log2, p={panel_correlation.pvalue:.3g}; "
              f"Spearman rho={spearmanr(panel_x, panel_y).statistic:.6f}; "
              f"{left_out} with no edited molecule left out")
        ax.plot(limits, limits, linestyle="--", color=REFERENCE_COLOR, linewidth=REFERENCE_WIDTH, zorder=0)
        ax.scatter(panel_x, panel_y, color=colors[library.enzyme], linewidths=0, clip_on=False,
                   **POINT_STYLE)
        ax.text(*CORRELATION_POSITION, f"r = {panel_correlation.statistic:.2f}",
                transform=ax.transAxes, va="top")
        # The panel names the points it draws. How many sequences were left out for carrying no edited
        # molecule is in the printed summary and the README, where there is room to say why.
        ax.text(*COUNT_POSITION, f"n = {len(panel)}", transform=ax.transAxes, va="bottom", ha="right")
        ax.set_title(title)
        ax.set(xlim=limits, ylim=limits)
        ax.set_xticks(ticks, labels=tick_labels)
        ax.set_yticks(ticks, labels=tick_labels)
        ax.set_box_aspect(1)
        if ax is not axes[0, 0]:
            ax.tick_params(left=False, labelleft=False)
            ax.spines["left"].set_visible(False)
    shared_x = figure.supxlabel(selection.x_label.replace("\\n", "\n"))
    place_shared_ylabel(figure, axes[0, 0], selection.y_label.replace("\\n", "\n"))
    centre_shared_xlabel(figure, shared_x)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    labels = [shared_x, figure.texts[-1] if figure.texts else shared_x]
    for ax in axes.flat:
        labels.extend([*ax.get_xticklabels(), *ax.get_yticklabels(), ax.title])
    for text in labels:
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            require(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width and box.y1 <= figure.bbox.height,
                    "A scatter label extends outside the page")
    for left, right in zip(axes.flat, list(axes.flat)[1:]):
        require(max(text.get_window_extent(renderer).x1 for text in left.get_xticklabels())
                <= min(text.get_window_extent(renderer).x0 for text in right.get_xticklabels()),
                "Neighbouring panels' tick labels overlap")
    save_figure(figure, selection.source)
    print("No in vivo placements are excluded; with/without-exclusion estimates are identical.")
