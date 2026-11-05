"""Per-variant agreement between direct and split recruitment of TadA to the boxB library.

Every axis carries one quantity built one way: a designed boxB sequence's mean A8/A2 editing as a ratio
to its matched wild type. Normalising to matched wild type removes the activity difference between
reactions at different enzyme concentrations and leaves the sequence determinants, which is what the
recruitment strategies are compared on. Every panel reads the same direct reaction on its y axis
against one split-recruitment reaction on its x axis, on one pair of equal axes.

Wild type is left out: it is the reference every ratio is taken against, so it would enter at 100% on
both axes by construction. A sequence with no edited molecule has no logarithm to place on these axes,
so it is left out and counted on the panel rather than given a floor value.

The panels are sized so that their r and n labels clear every point, which is measured below, and the
script fails rather than ship a figure where a label touches a marker.
"""

import math

import numpy as np
from scipy.stats import pearsonr, spearmanr

from analysis_helpers import (SAMPLE_KEYS, annotation, boxb_wildtype, designed_boxb_sequences,
                              require, selected_measurement, sequence_effects)
from plot_style import palettes, place_shared_ylabel, plt, save_pdf

# Panels 0.93 inch square, the size the cellular comparison draws, so the two can sit side by side.
FIGURE_SIZE = (2.56, 1.4)
# The marker the cellular comparison draws, so the two families read alike.
POINT_STYLE = dict(s=3.5, alpha=.7)
AXIS_PADDING = .3
# Labelled at 6.25, 25, and 100% of wild type wherever the axis ends, as the other variant scatters are.
TICK_POWERS = [-4, -2, 0]
REFERENCE_COLOR = "#999999"
REFERENCE_WIDTH = .4
CORRELATION_POSITION = (.04, .96)
COUNT_POSITION = (.96, .04)
# How far, in points, the edge of every marker must stand from the r and n labels.
LABEL_CLEARANCE = 1
AXIS_LABEL = "{enzyme} editing,\n% of wild type"
LOOP_REGION = "loop"

labels, _, colors = palettes()
selections = annotation("boxb_recruitment_comparison").sort_values("panel_order")
require(selections.notna().all().all() and selections.source.nunique() == 1 and selections.y_sample_id.nunique() == 1
        and selections.measurement.nunique() == 1 and not selections.x_sample_id.duplicated().any(),
        "A recruitment comparison is one figure reading one direct reaction against distinct split reactions")
source = selections.source.iloc[0]
direct = selections.y_sample_id.iloc[0]
# The reactions are read from the reaction annotation rather than from the two-reaction sets. This
# figure reads one direct reaction against several split ones, so the reactions it draws need not form
# a set, and which concentration of an editor it compares stays a choice of this selection rather than
# one inherited from the set another family draws. The palette names the enzymes, as everywhere else.
chosen = [direct, *selections.x_sample_id]
reactions = annotation("sample_info")
samples = reactions.loc[reactions.sample_id.isin(chosen), SAMPLE_KEYS[:4]].copy()
require(len(samples) == len(chosen), "Every compared reaction must be an annotated reaction")
samples["enzyme_label"] = samples.tada_type.map(labels)
samples["panel_order"] = range(len(samples))
require(samples.enzyme_label.notna().all(), "Every compared reaction needs an enzyme the palette names")
enzyme_of = samples.set_index("sample_id").tada_type


measurement = selected_measurement(selections.measurement.iloc[0])
effects, _ = sequence_effects(samples, measurement)
designed = designed_boxb_sequences()
effects = effects.merge(designed, on="hairpin", validate="many_to_one")
require(len(effects) == len(samples) * len(designed), "Every reaction must carry every designed sequence")
variants = effects.loc[effects.n_substitutions.gt(0)].pivot(
    index="hairpin", columns="sample_id", values="log2_pooled_ratio_to_matched_wt")
require(len(variants) == len(designed) - 1, "Incomplete designed variants")

# A window overlapping the loop holds boxB's own adenosines, so an edit there rewrites the genotype its
# molecule is counted as. The agreement is printed for each kind of window, since the pooled coefficient
# on the figure mixes them. A mutant's substitutions all lie in one window, so its first one places it.
regions, _, wildtype = boxb_wildtype()
loop = regions.loc[regions.region.eq(LOOP_REGION)].iloc[0]
windows = annotation("boxb_multimutant_heatmap_selection")
windows = windows.loc[windows.window_start.le(loop.end) & windows.window_end.ge(loop.start)]
first_change = variants.index.map(lambda hairpin: next(i for i, (a, b) in enumerate(zip(hairpin, wildtype), 1) if a != b))
loop_window = np.array([any(window.window_start <= p <= window.window_end for window in windows.itertuples()) for p in first_change])

panels = []
for selection in selections.itertuples():
    pairs = variants[[selection.x_sample_id, direct]]
    finite = np.isfinite(pairs).all(axis=1).to_numpy()
    panels.append((selection.x_sample_id, pairs.loc[finite], loop_window[finite], int((~finite).sum())))
values = np.concatenate([pairs.to_numpy().ravel() for _, pairs, _, _ in panels])
limits = (values.min() - AXIS_PADDING, values.max() + AXIS_PADDING)
tick_labels = [f"{100 * 2. ** power:.3g}" for power in TICK_POWERS]
require(limits[0] < min(TICK_POWERS) and limits[1] > max(TICK_POWERS), "Every labelled tick must sit on the axis")

figure, axes = plt.subplots(1, len(panels), figsize=FIGURE_SIZE, layout="constrained", sharex=True, sharey=True)
for ax, (sample_id, pairs, in_loop, dropped) in zip(axes, panels):
    x, y = pairs[sample_id].to_numpy(), pairs[direct].to_numpy()
    enzyme = labels[enzyme_of[sample_id]]
    correlation = pearsonr(x, y)
    # A coefficient tightens when effects compress, so the spread of both axes is printed beside it.
    print(f"{labels[enzyme_of[direct]]} {direct} (y) vs {enzyme} {sample_id} (x): n={len(x)}/{len(pairs) + dropped}, "
          f"{dropped} with no edited molecule left out; Pearson r={correlation.statistic:.6f} on log2, "
          f"p={correlation.pvalue:.3g}, r={pearsonr(100 * 2. ** x, 100 * 2. ** y).statistic:.6f} on percent; "
          f"Spearman rho={spearmanr(x, y).statistic:.6f}; SD x={x.std(ddof=1):.3f}, y={y.std(ddof=1):.3f} log2; "
          f"median x={100 * 2. ** np.median(x):.3f}%, y={100 * 2. ** np.median(y):.3f}% of WT")
    for name, mask in [("stem-window mutants", ~in_loop), ("loop-window mutants", in_loop)]:
        print(f"  {name}: n={mask.sum()}, Pearson r={pearsonr(x[mask], y[mask]).statistic:.6f} on log2; "
              f"SD x={x[mask].std(ddof=1):.3f}, y={y[mask].std(ddof=1):.3f} log2")
    ax.plot(limits, limits, linestyle="--", color=REFERENCE_COLOR, linewidth=REFERENCE_WIDTH, zorder=0)
    ax.scatter(x, y, color=colors[enzyme_of[sample_id]], linewidths=0, **POINT_STYLE)
    ax.text(*CORRELATION_POSITION, f"r = {correlation.statistic:.2f}", transform=ax.transAxes, va="top")
    # The panel names how many variants it draws. How many were left out, and why, is printed on every
    # run and stated in the README, which is where a caption can take it from.
    ax.text(*COUNT_POSITION, f"n = {len(x)}", transform=ax.transAxes, va="bottom", ha="right")
    ax.set(xlim=limits, ylim=limits, xlabel=AXIS_LABEL.format(enzyme=enzyme))
    ax.set_xticks(TICK_POWERS, labels=tick_labels)
    ax.set_yticks(TICK_POWERS, labels=tick_labels)
    ax.set_box_aspect(1)
    if ax is not axes[0]:
        ax.tick_params(left=False, labelleft=False)
        ax.spines["left"].set_visible(False)
place_shared_ylabel(figure, axes[0], AXIS_LABEL.format(enzyme=labels[enzyme_of[direct]]))

# Measure the drawn page: a label within reach of a marker, or text off the page, stops the figure here.
figure.canvas.draw()
renderer = figure.canvas.get_renderer()
marker_radius = math.sqrt(POINT_STYLE["s"]) / 2 * figure.dpi / 72
for ax, (sample_id, pairs, _, _) in zip(axes, panels):
    centres = ax.transData.transform(pairs[[sample_id, direct]].to_numpy())
    for text in ax.texts:
        box = text.get_window_extent(renderer)
        gap = np.hypot(centres[:, 0] - np.clip(centres[:, 0], box.x0, box.x1),
                       centres[:, 1] - np.clip(centres[:, 1], box.y0, box.y1)).min() - marker_radius
        require(gap >= LABEL_CLEARANCE * figure.dpi / 72, f"'{text.get_text()}' sits within {LABEL_CLEARANCE} pt of a marker in {sample_id}")
        print(f"  '{text.get_text()}' clears the nearest {sample_id} marker by {72 * gap / figure.dpi:.1f} pt")
    for text in [*ax.get_xticklabels(), *ax.get_yticklabels(), ax.xaxis.label]:
        box = text.get_window_extent(renderer)
        require(not text.get_text() or (box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width and box.y1 <= figure.bbox.height),
                "A label extends outside the page")
require(axes[0].xaxis.label.get_window_extent(renderer).x1 < axes[1].xaxis.label.get_window_extent(renderer).x0, "Neighbouring x labels overlap")
require(max(text.get_window_extent(renderer).x1 for text in axes[0].get_xticklabels())
        < min(text.get_window_extent(renderer).x0 for text in axes[1].get_xticklabels()), "Neighbouring panels' tick labels overlap")
box = axes[0].get_position()
print(f"Panels are {box.width * FIGURE_SIZE[0]:.2f} x {box.height * FIGURE_SIZE[1]:.2f} inch on a {FIGURE_SIZE[0]} x {FIGURE_SIZE[1]} inch page")
save_pdf(figure, source, png=True)
