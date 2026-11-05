"""Amino acid variant editing at two sites, pooling accepted reads across codons and barcodes.

Run from the repository root. Fractions use each site's own quality-callable denominator.
The two sites share reads and are not independent replicates. Derived values stay in memory.
"""

import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import gaussian_kde, pearsonr, spearmanr

from analysis_helpers import (FIGURE_MIN_LINKAGE_SUPPORT, MIN_CALLABLE_READS, POOLED_WT_LABEL,
                              comparison_sites, figure_sample, pooled_variant_sites, require)
from plot_style import plt, save_figure


OUTPUT = "amino_acid_site_correlation_A2_A15"

settings = {
    "size": (1.4122, 1.405), "point_size": 4, "point_alpha": .65, "color": "#AA4499",
    "log2_tick_step": 6, "log2_tick_offset": 0, "tick_format": r"$\mathdefault{2^{%d}}$",
    "reference_color": "#999999", "reference_width": .4, "annotation_position": (.04, .96),
    "axis_label": "% editing, {site}", "axis_label_pad": 3.24,
    "density_points": 256, "density_width_percent": 16, "density_gap_in": .04,
    "density_alpha": .6, "bandwidth_factor": .9, "iqr_to_sd": 1.34,
}

sample_id = figure_sample()
sites = comparison_sites()
pooled = pooled_variant_sites(sample_id, sites, POOLED_WT_LABEL,
                                FIGURE_MIN_LINKAGE_SUPPORT)
coverage = pooled.pivot(index="variant", columns="site_id", values="callable")
covered = coverage.ge(MIN_CALLABLE_READS).all(axis=1)
rates = 100 * pooled.pivot(index="variant", columns="site_id", values="fraction").loc[covered, sites]
require(len(rates) > 2 and np.isfinite(rates).all().all(), "Insufficient paired site coverage")
positive = rates.gt(0).all(axis=1)
points = np.log2(rates.loc[positive])
require(len(points) > 2, "Insufficient variants with positive editing at both sites")
x, y = points[sites[0]].to_numpy(), points[sites[1]].to_numpy()
require(np.std(x) > 0 and np.std(y) > 0, "Constant editing rates have no correlation")
correlation = pearsonr(x, y)
rank_correlation = spearmanr(x, y)
limits = (np.floor(min(x.min(), y.min())), np.ceil(max(x.max(), y.max())))
tick_step = settings["log2_tick_step"]
tick_start = np.ceil((limits[0] - settings["log2_tick_offset"]) / tick_step) * tick_step
ticks = np.arange(tick_start + settings["log2_tick_offset"], limits[1] + 1, tick_step)
tick_labels = [settings["tick_format"] % value for value in ticks]

figure, ax = plt.subplots(figsize=settings["size"], layout="constrained")
density_ax = make_axes_locatable(ax).append_axes(
    "right", size=f'{settings["density_width_percent"]}%', pad=settings["density_gap_in"], sharey=ax)
ax.plot(limits, limits, linestyle="--", color=settings["reference_color"],
        linewidth=settings["reference_width"], zorder=0)
ax.scatter(x, y, s=settings["point_size"], color=settings["color"],
           alpha=settings["point_alpha"], edgecolors="none", clip_on=False, zorder=3)
bandwidth = settings["bandwidth_factor"] * min(
    np.std(y, ddof=1), np.subtract(*np.percentile(y, [75, 25])) / settings["iqr_to_sd"]) * len(y) ** -.2
require(bandwidth > 0, "The marginal distribution has zero bandwidth")
density_y = np.linspace(*limits, settings["density_points"])
density = gaussian_kde(y, bw_method=bandwidth / np.std(y, ddof=1))(density_y)
density_ax.fill_betweenx(density_y, 0, density / density.max(), color=settings["color"],
                        alpha=settings["density_alpha"], linewidth=0)
density_ax.axvline(0, color=settings["reference_color"], linewidth=settings["reference_width"])
density_ax.set_xlim(0, 1)
density_ax.axis("off")
ax.text(*settings["annotation_position"], f"r = {correlation.statistic:.2f}",
        transform=ax.transAxes, ha="left", va="top")
ax.set(xlabel=settings["axis_label"].format(site=sites[0]), ylabel=settings["axis_label"].format(site=sites[1]),
       xlim=limits, ylim=limits, xticks=ticks, yticks=ticks)
ax.set_xticklabels(tick_labels)
ax.set_yticklabels(tick_labels)
ax.set_aspect("equal")
ax.xaxis.labelpad = ax.yaxis.labelpad = settings["axis_label_pad"]
figure.canvas.draw()
renderer = figure.canvas.get_renderer()
for label in [ax.xaxis.label, ax.yaxis.label, *ax.texts, *ax.get_xticklabels(), *ax.get_yticklabels()]:
    box = label.get_window_extent(renderer)
    require(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width
            and box.y1 <= figure.bbox.height, "A scatter label extends outside the page")
marker_radius = np.sqrt(settings["point_size"]) * figure.dpi / 144
point_pixels = ax.transData.transform(np.column_stack([x, y]))
for label in ax.texts:
    box = label.get_window_extent(renderer).padded(marker_radius)
    require(not any(box.contains(*point) for point in point_pixels), "A scatter label overlaps a variant")
save_figure(figure, OUTPUT)

print(f"{sample_id}: {len(rates)}/{len(coverage)} amino acid variants with callable reads at both sites; "
      f"{len(points)} positive at both sites plotted on log2 axes. "
      f"Pearson r on log2 rates={correlation.statistic:.8f}, p={correlation.pvalue:.6g}; "
      f"Spearman rho={rank_correlation.statistic:.8f}.")
print(f"Reference Pearson r on untransformed percentages for all covered variants, including zeros: "
      f"{pearsonr(*rates.to_numpy().T).statistic:.8f}.")
for site in sites:
    rows = pooled.loc[pooled.site_id.eq(site)]
    wt = rows.loc[rows.variant.eq(POOLED_WT_LABEL)]
    print(f"{site}: {rows.reads.sum():,} accepted reads, {rows.barcodes.sum():,} barcodes, "
          f"{rows.observed_codons.sum():,} observed codon variants; "
          f"{rows.edited.sum():,} edited / {rows['callable'].sum():,} callable reads. "
          f"Pooled synonymous WT: {wt.edited.item():,}/{wt['callable'].item():,} "
          f"= {100 * wt.fraction.item():.6f}%; zero-edited covered variants={(rates[site] == 0).sum()}.")
print(f"Undefined at one or both sites (omitted): {', '.join(coverage.index[~covered])}")
print(f"Zero editing at one or both sites ({(~positive).sum()} omitted from log axes): "
      f"{', '.join(rates.index[~positive])}")
print(f"All RNA reads from barcodes with >= {FIGURE_MIN_LINKAGE_SUPPORT} designed linkage reads "
      "pooled before division; no barcode or codon averaging, RNA coverage filter, "
      "wild-type normalization, deduplication, or pseudocount. Stop variants included; "
      "synonymous codons pooled into one WT point. Sites share reads and are not replicates.")
