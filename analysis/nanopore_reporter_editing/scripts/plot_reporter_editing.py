"""Draw manuscript Figure 1C from the reporter read counts."""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import fisher_exact

# Run from the repository root; scientific selections live in CSVs.
ANALYSIS = "analysis/nanopore_reporter_editing"
SOURCE = f"{ANALYSIS}/source_data/tada_lambdaN_replicates_tada.csv"
SAMPLES = f"{ANALYSIS}/annotations/figure_1c_samples.csv"
CLASSES = f"{ANALYSIS}/annotations/edit_classes.csv"
PALETTE_SOURCE = "analysis/boxb_in_vitro_sequencing/annotations/recorder_edit_summary_panels.csv"
PALETTE_COMPARISON = "recorder_edit_summary"
OUTPUT = f"{ANALYSIS}/figures/reporter_editing"
PLOT_DATA = "source_data/figure_1c.csv"
STATISTICS = f"{ANALYSIS}/tables/figure_1c_statistics.csv"
SAMPLE_COLUMN, DENOMINATOR_COLUMN = "Sample", "sequences without indels"
# The name used in the exported panel data and the statistics table.
LABEL_COLUMN = "sample_label"
TOTAL_COLUMN, INDEL_COLUMN = "Total sequences", "sequences with indels"
BOOTSTRAP_DRAWS, SEED, CONFIDENCE_LEVEL = 100_000, 20261004, .95
PERCENT_TOLERANCE = 1e-7
TEST_ALTERNATIVE = "two-sided"
STAR_THRESHOLDS, STAR_LABELS = (.001, .01, .05), ("***", "**", "*")
STAR_OFFSET_PT, FOLD_GAP_PT, ANNOTATION_TOP_PAD_PT, FOLD_DIGITS = 1, 1.5, 1, 1
FONT, FONT_SIZE, PNG_DPI = "Helvetica", 7, 96
FIGURE_SIZE, BAR_WIDTH_IN = (2.45, 2.00), .14
MARGINS_IN = {"left": .50, "right": .03, "top": .24, "bottom": .43}
PANEL_GAP_IN, PAIR_OFFSET, EDIT_GAP = .14, .20, 1.0
X_LIMITS = (-.5, 1.5)
Y_LABEL, ROW_LABELS = "% Edited RNA", ("boxB", "Edit count")
ROW_OFFSETS_PT, ROW_LABEL_GAP_PT, TITLE_PAD_PT = (-8, -20), 5, 5
AXIS_COLOR, AXIS_WIDTH, BAR_EDGE_WIDTH, ERROR_WIDTH = "#999999", .25, .2, .35
Y_LABEL_PAD_PT = 5
Y_LIMITS, Y_TICKS, CAP_SIZE = (0, 100), [0, 25, 50, 75, 100], 1.5
PADDING_INCHES = .02
PDF_METADATA = {"Creator": "Figure 1C count replot", "CreationDate": None, "ModDate": None}

font_manager.findfont(FONT, fallback_to_default=False)
plt.rcParams.update({
    "font.family": FONT, "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "axes.titlesize": FONT_SIZE, "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
    "legend.fontsize": FONT_SIZE, "legend.title_fontsize": FONT_SIZE,
    "axes.edgecolor": AXIS_COLOR, "axes.linewidth": AXIS_WIDTH,
    "xtick.major.width": AXIS_WIDTH, "ytick.major.width": AXIS_WIDTH,
    "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": .5,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "figure_1c",
})
source = pd.read_csv(SOURCE)
samples = pd.read_csv(SAMPLES, comment="#")
classes = pd.read_csv(CLASSES).sort_values("class_order")
palette_rows = pd.read_csv(PALETTE_SOURCE).set_index("name")
palette_row = palette_rows.loc[PALETTE_COMPARISON]
palettes = {column: dict(zip(palette_row.enzyme_labels.split(";"), palette_row[column].split(";")))
            for column in classes.palette_column.unique()}
assert all(set(samples.enzyme) <= set(palette) for palette in palettes.values()), "Missing Figure 2 enzyme colors"
assert source[SAMPLE_COLUMN].is_unique and samples[SAMPLE_COLUMN].is_unique
assert source[TOTAL_COLUMN].eq(source[INDEL_COLUMN] + source[DENOMINATOR_COLUMN]).all()
data = samples.merge(source, on=SAMPLE_COLUMN, how="left", validate="one_to_one", indicator=True)
assert data._merge.eq("both").all(), "A selected sample is missing from the source CSV"
assert not data.duplicated(["enzyme_order", "reporter_order"]).any()
assert data.groupby("enzyme_order").reporter_order.nunique().eq(samples.reporter_order.nunique()).all()
records = []
rng = np.random.default_rng(SEED)
tail = (1 - CONFIDENCE_LEVEL) / 2
for edit_class in classes.itertuples(index=False):
    for _, row in data.sort_values(["enzyme_order", "reporter_order"]).iterrows():
        k, n = int(row[edit_class.count_column]), int(row[DENOMINATOR_COLUMN])
        percent = 100 * k / n
        assert abs(percent - row[edit_class.percent_column]) <= PERCENT_TOLERANCE, "Selected percentage disagrees with counts"
        draws = rng.binomial(n, k / n, size=BOOTSTRAP_DRAWS)
        low, high = 100 * np.quantile(draws / n, [tail, 1 - tail])
        records.append({"sample": row[LABEL_COLUMN], "Sample": row[SAMPLE_COLUMN], "enzyme": row.enzyme,
                        "enzyme_order": row.enzyme_order, "reporter": row.reporter,
                        "reporter_order": row.reporter_order, "edit_class": edit_class.edit_class,
                        "class_order": edit_class.class_order, "edited_reads": k, "indel_free_reads": n,
                        "reported_percent": row[edit_class.percent_column], "percent": percent,
                        "ci_low": low, "ci_high": high, "bootstrap_draws": BOOTSTRAP_DRAWS, "seed": SEED})
plot_data = pd.DataFrame(records)
os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
os.makedirs(os.path.dirname(PLOT_DATA), exist_ok=True)
plot_data.drop(columns=[SAMPLE_COLUMN]).to_csv(PLOT_DATA, index=False)

enzymes = data[["enzyme_order", "enzyme"]].drop_duplicates().sort_values("enzyme_order")
display = plot_data.merge(samples[[SAMPLE_COLUMN, "plot_order", "boxb_label", "control"]],
                          on=SAMPLE_COLUMN, validate="many_to_one")

# Each comparison is conditional on observed read counts, with one sample per condition.
statistics = []
for (enzyme_order, class_order), group in display.groupby(["enzyme_order", "class_order"], sort=True):
    assert len(group) == 2 and group.control.sum() == 1, "Incomplete reporter/control comparison"
    present, absent = group.loc[~group.control].iloc[0], group.loc[group.control].iloc[0]
    assert absent.percent > 0, "Fold change needs a nonzero control fraction"
    table = [[int(row.edited_reads), int(row.indel_free_reads - row.edited_reads)] for row in (present, absent)]
    result = fisher_exact(table, alternative=TEST_ALTERNATIVE)
    statistics.append({"enzyme": present.enzyme, "enzyme_order": enzyme_order,
                       "edit_class": present.edit_class, "class_order": class_order,
                       "boxb_present_sample": present["sample"], "boxb_absent_sample": absent["sample"],
                       "present_edited": table[0][0], "present_other": table[0][1],
                       "absent_edited": table[1][0], "absent_other": table[1][1],
                       "fold_change_present_over_absent": present.percent / absent.percent,
                       "p_value": result.pvalue, "test": "Fisher exact", "alternative": TEST_ALTERNATIVE})
statistics = pd.DataFrame(statistics)
assert len(statistics) == len(enzymes) * len(classes)
statistics["significance"] = np.select([statistics.p_value.lt(cutoff) for cutoff in STAR_THRESHOLDS],
                                       STAR_LABELS, default="ns")
assert np.isfinite(statistics[["p_value", "fold_change_present_over_absent"]]).all().all()
statistics.to_csv(STATISTICS, index=False)

figure, axes = plt.subplots(1, len(enzymes), figsize=FIGURE_SIZE, sharey=True, squeeze=False)
panel_width_in = (FIGURE_SIZE[0] - MARGINS_IN["left"] - MARGINS_IN["right"]
                  - (len(enzymes) - 1) * PANEL_GAP_IN) / len(enzymes)
figure.subplots_adjust(left=MARGINS_IN["left"] / FIGURE_SIZE[0],
                       right=1 - MARGINS_IN["right"] / FIGURE_SIZE[0],
                       top=1 - MARGINS_IN["top"] / FIGURE_SIZE[1],
                       bottom=MARGINS_IN["bottom"] / FIGURE_SIZE[1],
                       wspace=PANEL_GAP_IN / panel_width_in)
bar_width = BAR_WIDTH_IN * (X_LIMITS[1] - X_LIMITS[0]) / panel_width_in
bar_labels, group_labels, row_labels, bar_positions = [], [], [], []
comparison_labels = []
for column, enzyme in enumerate(enzymes.itertuples(index=False)):
    ax = axes[0, column]
    for edit_class in classes.itertuples(index=False):
        center = edit_class.class_order * EDIT_GAP
        values = display.loc[display.enzyme_order.eq(enzyme.enzyme_order)
                             & display.class_order.eq(edit_class.class_order)].sort_values("plot_order")
        assert len(values) == 2 and values.control.sum() == 1, "Each edit class needs one reporter/control pair"
        for point in values.itertuples(index=False):
            x = center + (point.plot_order * 2 - 1) * PAIR_OFFSET
            color = edit_class.control_color if point.control else palettes[edit_class.palette_column][point.enzyme]
            ax.bar(x, point.percent, width=bar_width, color=color, edgecolor="black", linewidth=BAR_EDGE_WIDTH)
            ax.errorbar(x, point.percent, yerr=[[point.percent - point.ci_low], [point.ci_high - point.percent]],
                        color="black", fmt="none", linewidth=ERROR_WIDTH, capsize=CAP_SIZE, capthick=ERROR_WIDTH)
            label = ax.annotate(point.boxb_label, (x, 0), xycoords=("data", "axes fraction"),
                                xytext=(0, ROW_OFFSETS_PT[0]), textcoords="offset points", ha="center", va="top")
            bar_labels.append(label)
            bar_positions.append((ax, x, label))
        group_labels.append(ax.annotate(edit_class.edit_class, (center, 0), xycoords=("data", "axes fraction"),
                                        xytext=(0, ROW_OFFSETS_PT[1]), textcoords="offset points", ha="center", va="top"))
        result = statistics.loc[statistics.enzyme_order.eq(enzyme.enzyme_order)
                                & statistics.class_order.eq(edit_class.class_order)].iloc[0]
        top = values.ci_high.max()
        fold_offset = STAR_OFFSET_PT
        if result.significance != "ns":
            comparison_labels.append(ax.annotate(result.significance, (center, top),
                                                 xytext=(0, STAR_OFFSET_PT), textcoords="offset points",
                                                 ha="center", va="bottom"))
            fold_offset += FONT_SIZE + FOLD_GAP_PT
        comparison_labels.append(ax.annotate(f"{result.fold_change_present_over_absent:.{FOLD_DIGITS}f}x",
                                             (center, top), xytext=(0, fold_offset), textcoords="offset points",
                                             ha="center", va="bottom"))
    ax.set(xlim=X_LIMITS, ylim=Y_LIMITS, yticks=Y_TICKS, xticks=[])
    ax.set_title(enzyme.enzyme, pad=TITLE_PAD_PT)
    ax.tick_params(axis="y", color=AXIS_COLOR)
    if column == 0:
        ax.set_ylabel(Y_LABEL, labelpad=Y_LABEL_PAD_PT)
        for name, offset in zip(ROW_LABELS, ROW_OFFSETS_PT):
            row_labels.append(ax.annotate(name, (0, 0), xycoords="axes fraction",
                                          xytext=(-ROW_LABEL_GAP_PT, offset), textcoords="offset points",
                                          ha="right", va="top"))
    else:
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.spines["left"].set_visible(False)

# Fit one shared scale to the annotations, then check the actual rendered labels.
figure.canvas.draw()
renderer = figure.canvas.get_renderer()
upper = Y_LIMITS[1]
for label in comparison_labels:
    height_pt = label.axes.get_window_extent(renderer).height * 72 / figure.dpi
    reserve_pt = label.xyann[1] + label.get_window_extent(renderer).height * 72 / figure.dpi + ANNOTATION_TOP_PAD_PT
    assert reserve_pt < height_pt, "Statistical annotation does not fit the axes"
    upper = max(upper, label.xy[1] / (1 - reserve_pt / height_pt))
for ax in axes.flat:
    ax.set_ylim(Y_LIMITS[0], upper)
figure.canvas.draw()
renderer = figure.canvas.get_renderer()
labels = bar_labels + group_labels + row_labels + comparison_labels + [ax.title for ax in axes.flat]
for i, first in enumerate(labels):
    box = first.get_window_extent(renderer)
    assert figure.bbox.contains(box.x0, box.y0) and figure.bbox.contains(box.x1, box.y1), "Label outside figure"
    for second in labels[i + 1:]:
        assert not box.overlaps(second.get_window_extent(renderer)), "Overlapping condition labels"
for ax, x, label in bar_positions:
    box = label.get_window_extent(renderer)
    assert abs((box.x0 + box.x1) / 2 - ax.transData.transform((x, 0))[0]) < .1, "boxB label is not aligned to its bar"
assert sum(len(ax.patches) for ax in axes.flat) == len(plot_data), "Missing or duplicated plotted bar"
for extension in ["svg", "png", "pdf"]:
    metadata = PDF_METADATA if extension == "pdf" else {"Date": None} if extension == "svg" else None
    figure.savefig(f"{OUTPUT}.{extension}", dpi=PNG_DPI, bbox_inches="tight",
                   pad_inches=PADDING_INCHES, transparent=False, metadata=metadata)
plt.close(figure)
print(plot_data[["enzyme", "reporter", "edit_class", "edited_reads", "indel_free_reads", "percent", "ci_low", "ci_high"]].to_string(index=False))
print(statistics[["enzyme", "edit_class", "fold_change_present_over_absent", "p_value", "significance"]].to_string(index=False))
print(f"PDF/PNG size ratio: {os.path.getsize(f'{OUTPUT}.pdf') / os.path.getsize(f'{OUTPUT}.png'):.2f}; both retained")
