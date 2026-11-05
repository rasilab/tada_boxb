"""Show mean editing across the two sites beneath class-wise codon-variant box plots.

All classes retain distinct codons at each position; only barcodes of the same insert are pooled.
All peptide positions contribute. Each distinct DNA insert contributes equal weight.
A codon variant enters only where at least CLASS_MIN_BARCODES barcodes carry it at both sites: the long
low tail of this panel follows barcode support, not read depth, and one or two barcodes carry the
clone-to-clone spread rather than the variant's own editing.
Average the site percentages arithmetically after pooling counts within each DNA insert.
Axes use log2 spacing. Box statistics use raw mean percentages; zero means are omitted and counted.
"""

import numpy as np
import pandas as pd
from matplotlib.ticker import NullLocator
from scipy.stats import mannwhitneyu

from analysis_helpers import (CLASS_MIN_BARCODES, FIGURE_MIN_LINKAGE_SUPPORT, MIN_CALLABLE_READS,
                              POOLED_WT_LABEL, SIGNIFICANCE_LEVEL, comparison_sites, figure_sample,
                              pooled_variant_sites, require)
from plot_style import plt, save_figure


OUTPUT = "variant_class_distribution"
CLASSES = ["synonymous", "missense", "stop"]
TICK_LABELS = ["Synonymous", "Missense", "Nonsense"]
VARIANT_RESOLUTION = "codon"
BOX_WHISKER_IQR = 1.5
POINT_JITTER_SEED = 20261003
# Neighbouring classes are compared on the mean percentages drawn here, by rank, so the comparison
# needs no distributional assumption. Measured zero means enter the test and are only left off the
# logarithmic axis. Synonymous against nonsense follows from the two and is not drawn.
CLASS_COMPARISONS = [("synonymous", "missense"), ("missense", "stop")]
STAR_THRESHOLDS = [(.001, "***"), (.01, "**"), (SIGNIFICANCE_LEVEL, "*")]

settings = {
    "size": (1.2493, 1.61), "panel_aspect": 1,
    "point_colors": ["#888888", "#AA4499", "#D55E00"], "box_color": "#222222",
    "alpha": .45, "point_size": 3, "jitter_width": .22, "box_width": .5,
    "line_width": .7, "median_width": .9, "median_color": "#111111",
    "point_zorder": 1, "box_zorder": 4, "x_margin": .5,
    "percent_limits": (0, 100), "log2_tick_step": 3, "log2_tick_offset": 0,
    "tick_format": r"$\mathdefault{2^{%d}}$",
    "site_separator": "/", "y_label": "% editing, mean {sites}",
    "y_label_pad": 3, "y_label_y": .4,
    "tick_rotation": 45, "tick_alignment": "right", "tick_rotation_mode": "anchor",
    "bracket_octaves": 1.2, "bracket_headroom": 3.0, "bracket_line_width": .4, "bracket_inset": .14,
    "bracket_tick_octaves": .5, "bracket_color": "#222222", "bracket_star_pad": -1.5,
    "bracket_ns": "ns",
}

sample_id = figure_sample()
sites = comparison_sites()
classes, tick_labels = CLASSES, TICK_LABELS
rng = np.random.default_rng(POINT_JITTER_SEED)
require(len(classes) == len(tick_labels) == len(settings["point_colors"]), "Class labels and colors disagree")
pooled = pooled_variant_sites(sample_id, sites, POOLED_WT_LABEL,
                              FIGURE_MIN_LINKAGE_SUPPORT, variant_resolution=VARIANT_RESOLUTION)
require(pooled.designed_codons.eq(1).all() and pooled.observed_codons.le(1).all(),
        "Each distribution point must represent exactly one designed DNA insert")
coverage = pooled.pivot(index="variant", columns="site_id", values="callable")[sites]
rates = pooled.pivot(index="variant", columns="site_id", values="fraction")[sites]
support = pooled.pivot(index="variant", columns="site_id", values="barcodes")[sites]
covered = coverage.ge(MIN_CALLABLE_READS).all(axis=1) & support.ge(CLASS_MIN_BARCODES).all(axis=1)
rows = pooled[["variant", "variant_class"]].drop_duplicates().set_index("variant")
designed_counts = rows.variant_class.value_counts()
rows = rows.loc[covered & rows.variant_class.isin(classes)].copy()
rows["percent"] = 100 * rates.loc[rows.index].mean(axis=1, skipna=False)
require(rows.index.is_unique and len(rows) == covered.sum(), "Site averaging lost or duplicated a codon variant")
require(rows.percent.between(*settings["percent_limits"]).all(), "Editing rate outside the probability bounds")
positive = rows.loc[rows.percent.gt(0), "percent"]
require(len(positive) > 0, "No positive editing rate can be plotted on the log axis")
# Each class against each other class, on the covered variants the panel summarises, including the
# measured zeros the logarithmic axis cannot show.
comparisons = []
for first, second in CLASS_COMPARISONS:
    left, right = [rows.loc[rows.variant_class.eq(cls), "percent"] for cls in (first, second)]
    require(len(left) > 1 and len(right) > 1, f"Too few variants to compare {first} with {second}")
    result = mannwhitneyu(left, right, alternative="two-sided")
    comparisons.append(dict(first=first, second=second, n_first=len(left), n_second=len(right),
                            median_first=left.median(), median_second=right.median(),
                            statistic=result.statistic, p_value=result.pvalue))
comparisons = pd.DataFrame(comparisons)
comparisons["label"] = [next((symbol for limit, symbol in STAR_THRESHOLDS if p_value < limit),
                             settings["bracket_ns"]) for p_value in comparisons.p_value]

# The axis is ticked to the first step above the highest variant, and then ends above that tick to
# hold the comparison brackets.
data_limits = (np.floor(np.log2(positive.min())), np.ceil(np.log2(positive.max())))
tick_step = settings["log2_tick_step"]
tick_offset = settings["log2_tick_offset"]
first_tick = np.ceil((data_limits[0] - tick_offset) / tick_step) * tick_step + tick_offset
last_tick = first_tick + np.ceil((data_limits[1] - first_tick) / tick_step) * tick_step
tick_exponents = np.arange(first_tick, last_tick + 1, tick_step)
log_limits = (data_limits[0], last_tick + settings["bracket_headroom"])
limits = np.exp2(log_limits)
ticks = np.exp2(tick_exponents)
percent_labels = [settings["tick_format"] % value for value in tick_exponents]
class_positions = np.arange(1, len(classes) + 1)
figure, ax = plt.subplots(figsize=settings["size"], layout="constrained")
summaries = []
ax.set_box_aspect(settings["panel_aspect"])
ax.set_yscale("log", base=2)
ax.yaxis.set_minor_locator(NullLocator())
for position, cls, color in zip(class_positions, classes, settings["point_colors"]):
    observed = rows.loc[rows.variant_class.eq(cls), "percent"].to_numpy()
    values = observed[observed > 0]
    require(len(values) > 1, f"Insufficient positive variants for the mean-site {cls} box plot")
    jitter = rng.uniform(-settings["jitter_width"], settings["jitter_width"], len(values))
    ax.scatter(position + jitter, values, s=settings["point_size"], color=color,
               alpha=settings["alpha"], edgecolors="none", zorder=settings["point_zorder"])
    ax.boxplot([values], positions=[position], widths=settings["box_width"],
               whis=BOX_WHISKER_IQR, showfliers=False, patch_artist=True, manage_ticks=False,
               zorder=settings["box_zorder"],
               boxprops=dict(facecolor="none", edgecolor=settings["box_color"], linewidth=settings["line_width"]),
               whiskerprops=dict(color=settings["box_color"], linewidth=settings["line_width"]),
               capprops=dict(color=settings["box_color"], linewidth=settings["line_width"]),
               medianprops=dict(color=settings["median_color"], linewidth=settings["median_width"]))
    summaries.append(dict(variant_class=cls, designed_variants=designed_counts[cls],
                          covered_variants=len(observed), plotted_variants=len(values),
                          zero_omitted=(observed == 0).sum(), median_positive=np.median(values),
                          q1_positive=np.quantile(values, .25), q3_positive=np.quantile(values, .75),
                          minimum_positive=values.min(), maximum=values.max()))
# Brackets sit on one rung above the drawn data, each over the two classes it compares.
bracket_base = last_tick
for comparison in comparisons.itertuples():
    left, right = [class_positions[classes.index(cls)] for cls in (comparison.first, comparison.second)]
    # Two brackets meeting over the class they share would read as one, so each stops short of its ends.
    left, right = left + settings["bracket_inset"], right - settings["bracket_inset"]
    line = np.exp2(bracket_base + settings["bracket_octaves"])
    drop = np.exp2(bracket_base + settings["bracket_octaves"] - settings["bracket_tick_octaves"])
    ax.plot([left, left, right, right], [drop, line, line, drop],
            color=settings["bracket_color"], linewidth=settings["bracket_line_width"], clip_on=False)
    ax.annotate(comparison.label, ((left + right) / 2, line), xytext=(0, settings["bracket_star_pad"]),
                textcoords="offset points", ha="center", va="bottom", color=settings["bracket_color"])
ax.set(xlim=(class_positions.min() - settings["x_margin"], class_positions.max() + settings["x_margin"]),
       xticks=class_positions, xticklabels=tick_labels,
       ylim=limits, yticks=ticks, yticklabels=percent_labels,
       ylabel=settings["y_label"].format(sites=settings["site_separator"].join(sites)))
ax.yaxis.labelpad = settings["y_label_pad"]
ax.yaxis.label.set_y(settings["y_label_y"])
plt.setp(ax.get_xticklabels(), rotation=settings["tick_rotation"],
         ha=settings["tick_alignment"], rotation_mode=settings["tick_rotation_mode"])
figure.canvas.draw()
renderer = figure.canvas.get_renderer()
require(np.isclose(ax.bbox.height / ax.bbox.width, settings["panel_aspect"]),
        "The distribution panel does not have the requested aspect ratio")
texts = [ax.yaxis.label, *ax.get_xticklabels()]
texts.extend(label for tick, label in zip(ax.get_yticks(), ax.get_yticklabels()) if ax.get_ylim()[0] <= tick <= ax.get_ylim()[1])
boxes = [item.get_window_extent(renderer) for item in texts]
require(all(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width and box.y1 <= figure.bbox.height
            for box in boxes), "A distribution label extends outside the page")
rotated_labels = ax.get_xticklabels()
require(all(not left.overlaps(right)
            for index, left in enumerate(boxes) for next_index, right in enumerate(boxes[index + 1:], index + 1)
            if not (texts[index] in rotated_labels and texts[next_index] in rotated_labels)),
        "Distribution labels overlap")
# Rotated labels can have overlapping axis-aligned boxes while their text remains separate.
anchors = np.array([label.get_transform().transform(label.get_position()) for label in rotated_labels])
angle = np.deg2rad(settings["tick_rotation"])
normal_positions = np.sort(anchors @ np.array([-np.sin(angle), np.cos(angle)]))
require(np.diff(normal_positions).min() > max(label.get_fontsize() for label in rotated_labels) * figure.dpi / 72,
        "Rotated class labels are too close to read")
save_figure(figure, OUTPUT)
print(pd.DataFrame(summaries).to_string(index=False, float_format="%.6f"))
print("Class comparisons, two-sided Mann-Whitney on the drawn mean percentages including measured "
      f"zeros, uncorrected over the {len(comparisons)} pairs; *** is P < 0.001, ** P < 0.01, "
      f"* P < {SIGNIFICANCE_LEVEL:g}:")
print(comparisons.to_string(index=False, float_format="%.6g"))
print(f"Barcode floor: at least {CLASS_MIN_BARCODES} barcodes at both sites. "
      f"Covered codon variants fall from {int(coverage.ge(MIN_CALLABLE_READS).all(axis=1).sum())} to "
      f"{int(covered.sum())}; the discarded variants are those whose editing is carried by too few "
      "clones to separate from the spread wild-type clones themselves show.")
print(f"Variant resolution: {VARIANT_RESOLUTION} for every class; only barcodes of the same DNA insert are pooled. "
      "Distinct codons remain separate at every position; canonical WT appears once in the synonymous class. "
      f"Arithmetic mean of {'/'.join(sites)} percentages with paired coverage required; zero site rates retained in the mean. "
      "All peptide positions combined without averaging.")
print(f"Every positive-mean codon variant is drawn once, including box-plot outliers. Boxes use mean-percentage "
      f"quartiles and medians; whiskers extend to the last observation within {BOX_WHISKER_IQR:g} IQR. "
      "Horizontal jitter only separates points; log2 spacing preserves the untransformed editing rates. "
      "Zero means are omitted without a pseudocount or WT normalization.")
paired = coverage.ge(MIN_CALLABLE_READS).all(axis=1)
print("Undefined paired site rates (omitted): " + ", ".join(coverage.index[~paired]))
print(f"Below the {CLASS_MIN_BARCODES}-barcode floor (omitted): {int((paired & ~covered).sum())} codon variants, "
      f"carried by a median of {int(support.loc[paired & ~covered].min(axis=1).median())} barcodes.")
