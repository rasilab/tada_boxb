"""Accepted reads per codon variant, pooling barcodes for each designed DNA insert.

Synonymous codons remain separate; the exact wild-type DNA sequence appears once. Rank every
designed insert by decreasing reads. Omit zero-read variants from the logarithmic y-axis and
report their count without adding a pseudocount. Run from the repository root.
"""

import numpy as np
import pandas as pd
from matplotlib.ticker import NullLocator

from analysis_helpers import (COUNT_DIRECTORY, FIGURE_MIN_LINKAGE_SUPPORT, POOLED_WT_LABEL, figure_sample,
                              filter_linkage_support, inserts, require)
from plot_style import LABEL_GAP, plt, save_figure


OUTPUT = "variant_read_coverage"
Y_SCALE = "log"
Y_MIN, Y_MAX = 1, 1000000
Y_TICK_EXPONENT_STEP = 2
X_TICK_STEP = 500

settings = {"size": (1.4025, 1.33), "point_size": 1.5, "color": "#000000",
                          "wt_offset": (3, 0), "x_padding": .025}

sample_id = figure_sample()
designs = inserts()
designs["variant"] = designs.wt_codon + designs.codon_position.astype(str) + designs.mut_codon
designs.loc[designs.is_wild_type, "variant"] = POOLED_WT_LABEL
expected = len(designs)
require(designs.variant.is_unique, "Codon variant labels are not unique")
require(Y_SCALE == "log" and 0 < Y_MIN < Y_MAX, "Coverage requires positive logarithmic limits")
counts = pd.read_csv(f"{COUNT_DIRECTORY}/{sample_id}_reads.csv",
                     usecols=["barcode", "insert_num", "reads"])
require(len(counts) > 0 and counts.barcode.is_unique, "Empty or duplicated barcode counts")
require(counts.notna().all().all() and np.isfinite(counts.reads).all()
        and counts.reads.ge(0).all() and counts.reads.eq(np.floor(counts.reads)).all(),
        "Invalid accepted read counts")
require(counts.insert_num.isin(designs.row_num).all(), "Counts contain an unannotated insert")
counts = filter_linkage_support(counts, sample_id, FIGURE_MIN_LINKAGE_SUPPORT)
by_insert = counts.groupby("insert_num").reads.sum()
coverage = designs[["row_num", "variant"]].merge(
    by_insert, left_on="row_num", right_index=True, how="left", validate="one_to_one")
coverage["reads"] = coverage.reads.fillna(0).astype("int64")
require(len(coverage) == expected and coverage.reads.sum() == counts.reads.sum(),
        "Barcode pooling changed coverage or read totals")
coverage = coverage.sort_values(["reads", "variant"], ascending=[False, True]).reset_index(drop=True)
coverage["rank"] = np.arange(1, len(coverage) + 1)
wild_type = coverage.loc[coverage.variant.eq("WT")]
require(len(wild_type) == 1, "Wild type must appear exactly once")
observed = coverage.loc[coverage.reads.gt(0)]
uncovered = coverage.loc[coverage.reads.eq(0)]
require(len(observed) > 0, "No codon variant has reads")
require(observed.reads.between(Y_MIN, Y_MAX).all(),
        "A covered variant lies outside the declared y-axis limits")

figure, ax = plt.subplots(figsize=settings["size"], layout="constrained")
ax.scatter(observed["rank"], observed.reads, s=settings["point_size"],
           color=settings["color"], edgecolors="none", clip_on=False)
if wild_type.reads.item() > 0:
    ax.annotate("WT", (wild_type["rank"].item(), wild_type.reads.item()),
                xytext=settings["wt_offset"], textcoords="offset points", ha="left", va="center")
last_rank = observed["rank"].max()
ticks = [1, *range(X_TICK_STEP, last_rank + 1, X_TICK_STEP)]
ax.set(xlabel="Codon variant rank", ylabel="Reads per variant", xticks=ticks,
       xlim=(1 - settings["x_padding"] * last_rank, last_rank),
       yscale=Y_SCALE, ylim=(Y_MIN, Y_MAX))
powers = np.arange(np.log10(Y_MIN), np.log10(Y_MAX) + 1,
                   Y_TICK_EXPONENT_STEP).astype(int)
ax.set_yticks(10.0 ** powers, ["1" if power == 0 else f"10$^{{{power}}}$" for power in powers])
ax.yaxis.set_minor_locator(NullLocator())
ax.xaxis.labelpad = ax.yaxis.labelpad = LABEL_GAP * 72
figure.canvas.draw()
renderer = figure.canvas.get_renderer()
for label in [ax.xaxis.label, ax.yaxis.label, *ax.texts,
              *ax.get_xticklabels(), *ax.get_yticklabels()]:
    box = label.get_window_extent(renderer)
    require(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width
            and box.y1 <= figure.bbox.height, "A coverage label extends outside the page")
save_figure(figure, OUTPUT)

print(f"{sample_id}: {len(coverage)} codon variants including one exact WT; "
      f"{len(observed)} covered, {len(uncovered)} uncovered; {coverage.reads.sum():,} accepted reads.")
print(f"Reads per variant: minimum={coverage.reads.min():,}, "
      f"minimum positive={observed.reads.min():,}, median={coverage.reads.median():,.0f}, "
      f"maximum={coverage.reads.max():,}; WT={wild_type.reads.item():,}.")
print(f"Uncovered: {', '.join(uncovered.variant)}")
print("Counts include PCR duplicates; barcodes are pooled per DNA insert and synonymous codons remain separate.")
print(f"The logarithmic y-axis omits {len(uncovered)} zero-read variants; no pseudocount is added.")
