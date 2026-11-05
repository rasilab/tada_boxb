"""Show amino-acid editing distributions separately at each recorder site.

Gray boxes summarize equally weighted amino-acid rates, including zeros, after pooling codons
and barcodes. Whiskers use the percentiles declared below. These distributions describe variant spread,
not uncertainty in a pooled rate.
"""

import numpy as np
import pandas as pd
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from matplotlib.ticker import MultipleLocator

from analysis_helpers import (FIGURE_MIN_LINKAGE_SUPPORT, POOLED_WT_LABEL, figure_sample, layout,
                              measured_sites, pooled_variant_sites, require, reverse_complement, verify_sites)
from plot_style import FONT, FONT_SIZE, plt, save_figure


OUTPUT = "recorder_site_profile"
CLASSES = ["missense", "synonymous", "stop"]
VARIANT_RESOLUTION = "amino_acid"
WHISKER_LOW, WHISKER_HIGH = 5, 95

settings = {
    "size": (2.322, 1.24), "y_tick_step": 5, "y_limit_step": 1,
    "context_y": -.20, "position_y": -.38, "row_label_x": .48, "x_margin": .45,
    "y_label": "% editing", "y_label_pad": 3,
    "box_width": .42, "box_color": "#DDDDDD", "line_color": "#888888", "line_width": .5,
    "median_color": "#000000", "median_width": 1.2, "box_zorder": 1,
}

sample_id = figure_sample()
classes = CLASSES
primary_sites = measured_sites(primary_only=True)
require(primary_sites.feature.nunique() == 1, "The profile requires one recorder feature")
feature_name = primary_sites.feature.iloc[0]
require(0 <= WHISKER_LOW < 25 < 75 < WHISKER_HIGH <= 100, "Invalid profile whiskers")
verify_sites(sample_id)
feature = layout(sample_id)[feature_name]
require(feature.orientation == "antisense", "The counter requires an antisense recorder")
sequence = reverse_complement(feature.expected_sequence).replace("T", "U")
sites = measured_sites(feature=feature_name, primary_only=True).sort_values("rna_position").copy()
require(len(sites) > 0, "No recorder adenosines selected")
sites["derived_context"] = [sequence[position - 2:position + 1] for position in sites.rna_position]
require(sites.derived_context.eq(sites.context).all()
        and sites.derived_context.str[1].eq("A").all(), "Recorder contexts disagree with the sense sequence")
sites["position_index"] = np.arange(1, len(sites) + 1)
variants = pooled_variant_sites(sample_id, list(sites.site_id), POOLED_WT_LABEL,
                                FIGURE_MIN_LINKAGE_SUPPORT, variant_resolution=VARIANT_RESOLUTION)
require(set(classes) <= set(variants.variant_class), "Unknown distribution class")
variants = variants.loc[variants.variant_class.isin(classes)]
distributions = []
summaries = []
for site in sites.itertuples():
    site_variants = variants.loc[variants.site_id.eq(site.site_id)]
    rates = 100 * site_variants.loc[site_variants.callable.gt(0), "fraction"]
    require(len(rates) > 0 and rates.between(0, 100).all(), "Invalid amino-acid editing rates")
    low, q1, median, q3, high = np.percentile(
        rates, [WHISKER_LOW, 25, 50, 75, WHISKER_HIGH])
    distributions.append(dict(whislo=low, q1=q1, med=median, q3=q3, whishi=high))
    summaries.append(dict(site_id=site.site_id, covered=len(rates), zeros=int(rates.eq(0).sum()),
                          uncovered=len(site_variants) - len(rates), low=low, q1=q1, median=median, q3=q3, high=high,
                          minimum=rates.min(), maximum=rates.max(),
                          outside_whiskers=int(((rates < low) | (rates > high)).sum())))
figure, ax = plt.subplots(figsize=settings["size"], layout="constrained")
ax.bxp(distributions, positions=sites.position_index, widths=settings["box_width"],
       showfliers=False, patch_artist=True, manage_ticks=False, zorder=settings["box_zorder"],
       boxprops=dict(facecolor=settings["box_color"], edgecolor=settings["line_color"],
                     linewidth=settings["line_width"]),
       whiskerprops=dict(color=settings["line_color"], linewidth=settings["line_width"]),
       capprops=dict(color=settings["line_color"], linewidth=settings["line_width"]),
       medianprops=dict(color=settings["median_color"], linewidth=settings["median_width"]))
maximum = max(row["whishi"] for row in distributions)
maximum = np.ceil(maximum / settings["y_limit_step"]) * settings["y_limit_step"]
ax.set(ylim=(0, maximum), xlim=(1 - settings["x_margin"], len(sites) + settings["x_margin"]),
       ylabel=settings["y_label"], xticks=sites.position_index, xticklabels=[])
ax.yaxis.labelpad = settings["y_label_pad"]
ax.yaxis.set_major_locator(MultipleLocator(settings["y_tick_step"]))
ax.tick_params(axis="x", length=0)
labels = []
for site in sites.itertuples():
    context = HPacker(children=[TextArea(base, textprops=dict(fontfamily=FONT, fontsize=FONT_SIZE,
                       fontweight="bold" if index == 1 else "normal"))
                       for index, base in enumerate(site.derived_context)], align="baseline", pad=0, sep=0)
    artist = AnnotationBbox(context, (site.position_index, settings["context_y"]),
                           xycoords=ax.get_xaxis_transform(), box_alignment=(.5, .5),
                           frameon=False, pad=0, annotation_clip=False)
    ax.add_artist(artist)
    labels.append(artist)
    labels.append(ax.text(site.position_index, settings["position_y"], str(site.rna_position),
                          ha="center", va="center", transform=ax.get_xaxis_transform()))
for y, label in [(settings["context_y"], "A context"), (settings["position_y"], "A position")]:
    labels.append(ax.text(settings["row_label_x"], y, label, ha="right", va="center",
                          transform=ax.get_xaxis_transform()))
figure.canvas.draw()
renderer = figure.canvas.get_renderer()
boxes = [label.get_window_extent(renderer) for label in [*labels, ax.yaxis.label]]
require(all(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width and box.y1 <= figure.bbox.height
            for box in boxes), "A recorder profile label extends outside the page")
require(all(not left.overlaps(right) for index, left in enumerate(boxes) for right in boxes[index + 1:]),
        "Recorder profile labels overlap")
save_figure(figure, OUTPUT)
print(f"{sample_id}: linkage support >= {FIGURE_MIN_LINKAGE_SUPPORT}; "
      "codon and barcode counts pooled per amino-acid variant before division by site-callable reads.")
print(f"Amino-acid distributions include {';'.join(CLASSES)} and zero-edit rates, with equal weight "
      f"per covered variant and one pooled synonymous {POOLED_WT_LABEL} reference. "
      f"Boxes = 25th–75th percentiles with a median line; whiskers = "
      f"{WHISKER_LOW:g}th–{WHISKER_HIGH:g}th percentiles. "
      "Rates outside the whiskers are summarized below but not drawn as outliers.")
print(pd.DataFrame(summaries).to_string(index=False, float_format="%.6f"))
print("Site coordinates and sense contexts verified. Boxes show observed variation among amino-acid "
      "variants, including counting noise; they are not confidence intervals or biological-replicate error bars.")
