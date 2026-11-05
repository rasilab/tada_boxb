"""Equal-height amino-acid logos of arithmetic mean site editing at each peptide position.

Pool codon/barcode counts per amino acid before calculating each site's rate. Average the
declared sites, including zeros, then divide by the summed amino-acid rates within a position.
WT is a shared measured reference; missing rates remain undefined and are never set to zero.
A letter is outlined where its own substitution differs from the pooled wild type, and a position
is starred where the mean of its substitutions does. Both marks require the two sites to agree and
count off the least favourable of their two false discovery rates, as the position profiles do; an
unmarked letter is as often a variant too few barcodes carry as one that matches wild type.
"""

import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import Patch, PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

from analysis_helpers import (FIGURE_MIN_LINKAGE_SUPPORT, MIN_CALLABLE_READS, POOLED_WT_LABEL,
                              SIGNIFICANCE_LEVEL, adjust_bh, annotation, comparison_sites, figure_sample,
                              inserts, pooled_variant_sites, position_mean_test, require,
                              wild_type_comparisons)
from plot_style import FONT, plt, save_figure

OUTPUT = "amino_acid_editing_logo"
CLASSES = ["missense"]
RESIDUE_ORDER = "AVLIMFWGPSTCYNQKRHDE"
GROUP_LABELS = {"hydrophobic": "Hydrophobic", "polar": "Polar", "basic": "Basic", "acidic": "Acidic"}

# Retain the 3.50-inch axes width after allowing for the legend and labels.
settings = {
    "size": (4.87193, 1.596), "letter_width": .82, "letter_weight": "bold",
    "group_colors": {"hydrophobic": "#888888", "polar": "#44AA99", "basic": "#4477AA", "acidic": "#D55E00"},
    "wt_color": "#222222", "wt_label": "WT residue",
    "variant_alpha": .60, "wt_alpha": 1,
    "x_margin": .5, "position_label_step": 2, "y_limits": (0, 1), "y_ticks": (0, .5, 1),
    "site_separator": "/",
    "x_label": "Position along λN",
    "y_label": "Relative editing\nmean {sites}", "y_label_pad": 3, "y_label_y": .5,
    "legend_handle_length": 1,
    "legend_label_spacing": .5, "legend_handletextpad": .4,
    "marked_edge_color": "#000000", "marked_edge_width": .35,
    "star_thresholds": [(.001, "***"), (.01, "**"), (.05, "*")], "star_y": 1.02,
}


def logo_letter(ax, letter, x, bottom, width, height, color, alpha=1, outlined=False):
    """Draw a vector Helvetica glyph with height encoding a normalized editing rate."""
    path = TextPath((0, 0), letter, size=1,
                    prop=font_manager.FontProperties(family=FONT, weight=settings["letter_weight"]))
    box = path.get_extents()
    require(box.width > 0 and box.height > 0 and height > 0, f"Invalid logo glyph for {letter}")
    glyph_width = width * min(1, box.width / box.height)
    transform = Affine2D().translate(-box.x0, -box.y0).scale(glyph_width / box.width, height / box.height).translate(x - glyph_width / 2, bottom)
    glyph = PathPatch(transform.transform_path(path), facecolor=color, alpha=alpha,
                      edgecolor=settings["marked_edge_color"] if outlined else "none",
                      linewidth=settings["marked_edge_width"] if outlined else 0)
    ax.add_patch(glyph)
    return glyph


sample_id = figure_sample()
sites = comparison_sites()
classes = CLASSES
designs = inserts()
order = annotation("amino_acid_groups")
order["display_order"] = order.amino_acid.map({residue: index for index, residue in enumerate(RESIDUE_ORDER)})
order["label"] = order.amino_acid
order["group_label"] = order.group.map(GROUP_LABELS)
order = order.sort_values("display_order")
residues = designs[["codon_position", "wt_aa"]].drop_duplicates().sort_values("codon_position")
positions = residues.codon_position.to_numpy()
amino_acids = designs.loc[designs.variant_class.isin([*classes, "synonymous"]), "mut_aa"].unique()
order_for_logo = order.loc[order.amino_acid.isin(amino_acids)]
require(order_for_logo.amino_acid.is_unique and set(order_for_logo.amino_acid) == set(amino_acids)
        and order_for_logo.group.isin(settings["group_colors"]).all()
        and order_for_logo[["display_order", "group_label"]].notna().all().all(),
        "Missing amino-acid group or display setting")

pooled = pooled_variant_sites(sample_id, sites, POOLED_WT_LABEL,
                                FIGURE_MIN_LINKAGE_SUPPORT)
# One test per amino-acid variant at each site, against the pooled synonymous wild type at the
# dispersion its own barcodes show. A mark needs both sites, as it does on the position profiles.
comparisons, reference = wild_type_comparisons(sample_id, sites, POOLED_WT_LABEL,
                                               FIGURE_MIN_LINKAGE_SUPPORT)
variant_rates = comparisons.pivot(index="variant", columns="site_id", values="q_value_barcode")[sites]
marked = set(variant_rates.index[variant_rates.lt(SIGNIFICANCE_LEVEL).all(axis=1)])
coverage = pooled.pivot(index="variant", columns="site_id", values="callable")[sites]
rates = pooled.pivot(index="variant", columns="site_id", values="fraction")[sites]
covered = coverage.ge(MIN_CALLABLE_READS).all(axis=1)
require(covered.loc[POOLED_WT_LABEL], "The synonymous WT reference lacks paired site coverage")
cells = residues.merge(order_for_logo, how="cross").rename(columns={"amino_acid": "mut_aa"})
cells["is_reference"] = cells.mut_aa.eq(cells.wt_aa)
cells["variant"] = cells.wt_aa + cells.codon_position.astype(str) + cells.mut_aa
cells.loc[cells.is_reference, "variant"] = POOLED_WT_LABEL
designed_mutants = designs.loc[designs.variant_class.isin(classes)]
expected_variants = set(designed_mutants.wt_aa + designed_mutants.codon_position.astype(str) + designed_mutants.mut_aa)
require(set(cells.loc[~cells.is_reference, "variant"]) == expected_variants
        and cells.variant.isin(rates.index).all(), "Logo grid does not match the designed amino-acid variants")
cells = cells.merge(rates.reset_index(), on="variant", how="left", validate="many_to_one")
cells["covered"] = cells.variant.map(covered)
cells["marked"] = cells.variant.isin(marked) & ~cells.is_reference
cells["mean_editing"] = cells[sites].mean(axis=1, skipna=False).where(cells.covered)
position_totals = cells.groupby("codon_position").mean_editing.sum(min_count=1)
require(position_totals.gt(0).all() and np.isfinite(position_totals).all(),
        "A position has no positive observed editing to normalize")
cells["height"] = cells.mean_editing / cells.codon_position.map(position_totals)
require(np.allclose(cells.groupby("codon_position").height.sum(), 1)
        and cells.loc[~cells.covered, "height"].isna().all()
        and cells.loc[cells.mean_editing.eq(0), "height"].eq(0).all(),
        "Logo normalization lost unit height or changed missing/zero values")
# The same test at the level the stacks are read: the equal-weight mean over a position's covered
# substitutions against the wild-type rate, at each site, corrected over the 22 positions.
substitutions = cells.loc[~cells.is_reference & cells.covered, ["codon_position", "variant"]].copy()
position_tests = []
for site in sites:
    site_rate = comparisons.loc[comparisons.site_id.eq(site)].set_index("variant").rate
    row = reference.set_index("site_id").loc[site]
    for position, group in substitutions.groupby("codon_position"):
        _, statistic, p_value = position_mean_test(group.variant.map(site_rate), row.reference_rate,
                                                   row.reference_error)
        position_tests.append(dict(codon_position=position, site_id=site, statistic=statistic,
                                   substitutions=int(group.variant.map(site_rate).notna().sum()),
                                   p_value=p_value))
position_tests = pd.DataFrame(position_tests)
position_tests["q_value"] = np.nan
for site in sites:
    at_site = position_tests.site_id.eq(site)
    position_tests.loc[at_site, "q_value"] = adjust_bh(position_tests.loc[at_site, "p_value"]).to_numpy()
position_rates = position_tests.pivot(index="codon_position", columns="site_id", values="q_value")[sites]
require(len(position_rates) == len(positions) and position_rates.notna().all().all(),
        "A position lacks a test against wild type")
starred = position_rates.max(axis=1).loc[position_rates.lt(SIGNIFICANCE_LEVEL).all(axis=1)]
stars = {position: next(symbol for limit, symbol in settings["star_thresholds"] if q_value < limit)
         for position, q_value in starred.items()}

missing = cells.loc[~cells.covered]
zeros = cells.loc[cells.covered & cells.mean_editing.eq(0)]
positive = cells.loc[cells.height.gt(0)]

figure, ax = plt.subplots(figsize=settings["size"], layout="constrained")
glyphs = []
for position in positions:
    stack = positive.loc[positive.codon_position.eq(position)].sort_values(["height", "display_order"])
    bottom = 0
    for row in stack.itertuples():
        color = settings["wt_color"] if row.is_reference else settings["group_colors"][row.group]
        alpha = settings["wt_alpha"] if row.is_reference else settings["variant_alpha"]
        glyph = logo_letter(ax, row.label, position, bottom, settings["letter_width"], row.height,
                            color, alpha, outlined=row.marked)
        bounds = glyph.get_path().get_extents()
        require(np.isclose(bounds.height, row.height) and np.isclose(bounds.y0, bottom),
                "A rendered letter disagrees with its normalized editing rate")
        glyphs.append(glyph)
        bottom += row.height
    require(np.isclose(bottom, 1), "A rendered logo stack does not reach unit height")
require(len(glyphs) == len(positive), "Positive amino-acid variants were lost or duplicated in rendering")
for position, symbol in stars.items():
    ax.text(position, settings["star_y"], symbol, ha="center", va="bottom",
            color=settings["marked_edge_color"], transform=ax.get_xaxis_transform())
groups = order_for_logo[["group", "group_label"]].drop_duplicates()
handles = [Patch(facecolor=settings["group_colors"][row.group], label=row.group_label,
                 alpha=settings["variant_alpha"]) for row in groups.itertuples()]
handles.append(Patch(facecolor=settings["wt_color"], label=settings["wt_label"], alpha=settings["wt_alpha"]))
legend = figure.legend(handles=handles, loc="outside right center", frameon=False,
                       handlelength=settings["legend_handle_length"], labelspacing=settings["legend_label_spacing"],
                       handletextpad=settings["legend_handletextpad"])
ax.set(xlim=(positions.min() - settings["x_margin"], positions.max() + settings["x_margin"]),
       ylim=settings["y_limits"], yticks=settings["y_ticks"],
       ylabel=settings["y_label"].format(sites=settings["site_separator"].join(sites)), xlabel=settings["x_label"])
ax.set_xticks(positions[::settings["position_label_step"]])
ax.yaxis.labelpad = settings["y_label_pad"]
ax.yaxis.label.set_y(settings["y_label_y"])
figure.canvas.draw()
# Keep the measured panel width fixed in both PDF and PNG exports.
figure.set_layout_engine("none")
renderer = figure.canvas.get_renderer()
texts = [ax.xaxis.label, ax.yaxis.label, legend, *ax.get_xticklabels(), *ax.get_yticklabels()]
boxes = [item.get_window_extent(renderer) for item in texts]
require(all(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width and box.y1 <= figure.bbox.height
            for box in boxes), "A logo label extends outside the page")
require(all(not left.overlaps(right) for index, left in enumerate(boxes) for right in boxes[index + 1:]),
        "Logo labels overlap")
save_figure(figure, OUTPUT)

report = cells.groupby("codon_position").agg(covered_amino_acids=("covered", "sum"),
                                            total_editing=("mean_editing", "sum"), stack_height=("height", "sum"))
report["positive_letters"] = positive.groupby("codon_position").size()
report["zero_amino_acids"] = zeros.groupby("codon_position").size().reindex(positions, fill_value=0)
report["missing_amino_acids"] = missing.groupby("codon_position").size().reindex(positions, fill_value=0)
top = positive.sort_values(["codon_position", "height", "display_order"]).groupby("codon_position").tail(1).set_index("codon_position")
report["highest_editing_amino_acid"] = top.mut_aa
report["highest_mean_editing_percent"] = 100 * top.mean_editing
print(report.to_string(float_format="%.6f"))
print(f"Mean sites: {';'.join(sites)}; arithmetic average after pooling codons/barcodes per amino acid. "
      "Measured zero rates remain in the site mean. No stops, pseudocounts, or information-content transform.")
print("Letter height = mean editing / summed mean editing among covered amino acids at that position. "
      "Every stack totals one; absolute editing across positions cannot be compared from stack heights. "
      "Letters increase in editing from bottom to top, with annotation order breaking exact ties.")
print(f"Shared synonymous WT reference ({POOLED_WT_LABEL}): "
      f"{(100 * rates.loc[POOLED_WT_LABEL]).to_dict()}; "
      "repeated only as a reference for the WT residue at each position, not independent measurements.")
print(f"Covered missense variants: {cells.loc[cells.covered & ~cells.is_reference, 'variant'].nunique()}; "
      f"WT reference placements: {cells.is_reference.sum()}; positive letters: {len(positive)}.")
print("Missing amino-acid measurements: " + ", ".join(missing.variant))
print("Measured zero means (zero-height letters): " + ", ".join(zeros.variant))
drawn = positive.loc[positive.marked, "variant"].nunique()
print(f"Outlined letters, Benjamini-Hochberg q < {SIGNIFICANCE_LEVEL:g} at both sites against the "
      "pooled wild type ("
      + "; ".join(f"{row.site_id}={100 * row.reference_rate:.2f}% at dispersion {row.dispersion:.4f}"
                  for row in reference.itertuples())
      + f" over {reference.reference_barcodes.iloc[0]} barcodes): {drawn} of "
      f"{positive.loc[~positive.is_reference, 'variant'].nunique()} drawn substitutions. "
      "The test counts the barcodes carrying a variant, not its reads, so an unmarked letter is as "
      "often a variant too few barcodes carry as one that matches wild type.")
print("Marked substitutions: " + ", ".join(sorted(positive.loc[positive.marked, "variant"].unique())))
print(position_tests.to_string(index=False, float_format="%.6f"))
print("Starred positions, equal-weight mean over their substitutions, both sites agreeing, at the "
      "least favourable of their two rates: "
      + ", ".join(f"{position} {symbol} (q={starred[position]:.4f})" for position, symbol in stars.items())
      + ". *** is q < 0.001, ** q < 0.01, and * q < 0.05.")
