"""Plot editing sites separately, weighting amino acid variants equally at each position.

Run from the repository root. Codon/barcode counts are pooled before calculating each site rate.
Zero editing contributes to arithmetic means; undefined site rates never become zero.
Each position is tested against the pooled synonymous wild type at each site separately, and is
starred where both sites agree, at the least favourable of their two false discovery rates.
"""

import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator

from analysis_helpers import (FIGURE_MIN_LINKAGE_SUPPORT, MIN_CALLABLE_READS, POOLED_WT_LABEL,
                              SIGNIFICANCE_LEVEL, adjust_bh, clustered_rate_test, comparison_sites,
                              figure_sample, inserts, pooled_variant_sites, position_mean_test, require,
                              variant_barcode_sites, wild_type_reference)
from plot_style import plt, save_figure


PROFILES = [
    ("missense", "position_editing_profile", "Mean % editing\nmissense variants", (3.12925, 1.4086), 15),
]

settings = {
    "colors": ["#AA4499", "#117733"], "markers": ["o", "s"],
    "marker_size": 4, "line_width": .5,
    "y_tick_step": 5, "x_margin": .5, "reference_width": .5, "y_label_pad": 3,
    "reference_alpha": .85, "reference_style": ":", "reference_zorder": 0, "data_zorder": 3,
    "legend_handle_length": 1.1, "legend_label_spacing": .4, "legend_handletextpad": .35,
    "legend_borderpad": 0, "legend_borderaxespad": 0,
    "legend_location": "outside right center",
    "inset_legend_location": "lower right", "inset_legend_anchors": {"missense": (.985, .04)},
    "position_label_step": {"missense": 2},
    "position_label_format": {"missense": "{wt}\n{position}"},
    "legend_columns": 2, "legend_column_spacing": .5, "legend_site_suffix": ",",
    "star_thresholds": [(.001, "***"), (.01, "**"), (.05, "*")], "star_offset": -7,
    # Positions 7 and 8 share the same lower point, so their stars would sit on one line and read as
    # one mark; position 8 drops as far below as the panel floor allows.
    "star_extra_offset": {"missense": {8: -3}},
    "reference_label": "WT", "y_label_y": .5, "x_label": "Position along λN",
}

sample_id = figure_sample()
sites = comparison_sites()
designs = inserts()
# The barcodes of one variant disagree far more than their read counts imply, so the wild-type
# barcodes set both the reference rate every test uses and the dispersion it is tested at.
barcode_counts, _, _ = variant_barcode_sites(sample_id, sites, POOLED_WT_LABEL, FIGURE_MIN_LINKAGE_SUPPORT)
wild_type = wild_type_reference(barcode_counts, POOLED_WT_LABEL).set_index("site_id")

for variant_class, output, y_label, size, maximum in PROFILES:
    classes = [variant_class]
    residues = designs[["codon_position", "wt_aa"]].drop_duplicates().sort_values("codon_position")
    positions = residues.codon_position.to_numpy()

    pooled = pooled_variant_sites(sample_id, sites, POOLED_WT_LABEL,
                                    FIGURE_MIN_LINKAGE_SUPPORT)
    coverage = pooled.pivot(index="variant", columns="site_id", values="callable")[sites]
    rates = pooled.pivot(index="variant", columns="site_id", values="fraction")[sites]
    covered = coverage.ge(MIN_CALLABLE_READS).all(axis=1)
    # Both sites use the same covered variants; no averaging across sites.
    wt_percent = 100 * rates.loc[POOLED_WT_LABEL]
    require(covered.loc[POOLED_WT_LABEL] and np.isfinite(wt_percent).all(),
            "The pooled synonymous WT lacks paired site coverage")
    cells = designs.loc[designs.variant_class.isin(classes),
                        ["codon_position", "wt_aa", "mut_aa", "variant_class"]].drop_duplicates().copy()
    cells["variant"] = cells.wt_aa + cells.codon_position.astype(str) + cells.mut_aa
    require(cells.variant.is_unique and cells.variant.isin(rates.index).all(),
            "Position-specific amino acid variants are repeated or missing from pooled counts")
    cells = cells.merge(pooled[["variant", "site_id", "fraction"]], on="variant", validate="one_to_many")
    cells["percent"] = 100 * cells.fraction.where(cells.variant.map(covered))
    profile = cells.groupby(["codon_position", "site_id"]).agg(
        percent=("percent", "mean"), covered_variants=("percent", "count"),
        designed_variants=("percent", "size"), zero_variants=("percent", lambda values: values.eq(0).sum())).reset_index()
    profile = residues.merge(profile, on="codon_position", how="left", validate="one_to_many")
    require(profile.percent.notna().all() and profile.percent.between(0, 100).all(),
            "A peptide position lacks a valid mean editing rate")
    require(len(profile) == len(positions) * len(sites)
            and profile.covered_variants.sum() == cells.percent.notna().sum(),
            "Position averaging lost or duplicated an amino acid variant")

    # A missense position is tested on the substitutions the panel averages, so its spread carries
    # the barcode noise inside each of them; a nonsense position pools one group, so its own
    # barcodes are the units. False discovery rates are taken over the positions of one site, since
    # the two sites are measured on the same molecules and are not independent tests.
    tests = []
    for site in sites:
        row = wild_type.loc[site]
        for position, group in cells.loc[cells.site_id.eq(site)].groupby("codon_position"):
            if variant_class == "missense":
                _, statistic, p_value = position_mean_test(group.percent / 100, row.reference_rate,
                                                           row.reference_error)
                units = int(group.percent.notna().sum())
            else:
                carriers = barcode_counts.loc[barcode_counts.variant.isin(group.variant)
                                              & barcode_counts.site_id.eq(site)
                                              & barcode_counts["callable"].gt(0)]
                _, statistic, p_value = clustered_rate_test(carriers.edited, carriers["callable"],
                                                            row.reference_rate, row.reference_error,
                                                            row.dispersion)
                units = int(len(carriers))
            tests.append(dict(codon_position=position, site_id=site, tested_units=units,
                              statistic=statistic, p_value=p_value))
    tests = pd.DataFrame(tests)
    tests["q_value"] = np.nan
    for site in sites:
        at_site = tests.site_id.eq(site)
        tests.loc[at_site, "q_value"] = adjust_bh(tests.loc[at_site, "p_value"]).to_numpy()
    tests["significant"] = tests.q_value.lt(SIGNIFICANCE_LEVEL)
    profile = profile.merge(tests, on=["codon_position", "site_id"], validate="one_to_one")
    require(profile.p_value.notna().all(), "A position lacks a test against wild type")
    # A position is marked only where both sites agree, and its stars count off the least favourable
    # of the two rates, so the mark never claims more than the weaker measurement supports.
    rates = tests.pivot(index="codon_position", columns="site_id", values="q_value")[sites]
    marks = rates.max(axis=1).loc[rates.lt(SIGNIFICANCE_LEVEL).all(axis=1)]
    stars = {position: next(symbol for limit, symbol in settings["star_thresholds"] if q_value < limit)
             for position, q_value in marks.items()}

    figure, ax = plt.subplots(figsize=size, layout="constrained")
    site_handles = []
    reference_handles = []
    for index, site in enumerate(sites):
        rows = profile.loc[profile.site_id.eq(site)].sort_values("codon_position")
        color, marker = settings["colors"][index], settings["markers"][index]
        line, = ax.plot(rows.codon_position, rows.percent, marker=marker, color=color,
                       markersize=settings["marker_size"], linewidth=settings["line_width"],
                       markeredgewidth=0, clip_on=False, label=site, zorder=settings["data_zorder"])
        reference = ax.axhline(wt_percent[site], color=color, linewidth=settings["reference_width"],
                              linestyle=settings["reference_style"], alpha=settings["reference_alpha"],
                              zorder=settings["reference_zorder"])
        site_handles.append(line)
        reference_handles.append(reference)
    inset_legend = variant_class in settings["inset_legend_anchors"]
    legend_options = dict(loc=settings["legend_location"])
    if inset_legend:
        legend_options = dict(loc=settings["inset_legend_location"],
                              bbox_to_anchor=settings["inset_legend_anchors"][variant_class],
                              bbox_transform=ax.transAxes)
    # The stars sit under the lower of the two sites, in the space a depressed position leaves.
    floor = profile.groupby("codon_position").percent.min()
    for position, symbol in stars.items():
        offset = settings["star_offset"] + settings["star_extra_offset"][variant_class].get(position, 0)
        ax.annotate(symbol, (position, floor[position]), xytext=(0, offset),
                    textcoords="offset points", ha="center", va="top", zorder=settings["data_zorder"])
    legend = figure.legend(handles=site_handles + reference_handles,
                       labels=[f"{site}{settings['legend_site_suffix']}" for site in sites]
                              + [settings["reference_label"]] * len(sites),
                       **legend_options,
                       ncols=settings["legend_columns"], columnspacing=settings["legend_column_spacing"],
                       frameon=False, handlelength=settings["legend_handle_length"],
                       labelspacing=settings["legend_label_spacing"],
                       handletextpad=settings["legend_handletextpad"], borderpad=settings["legend_borderpad"],
                       borderaxespad=settings["legend_borderaxespad"])
    require(max(profile.percent.max(), wt_percent.max()) <= maximum, "The y-axis clips a position or WT mean")
    ax.set(ylim=(0, maximum), xlim=(positions.min() - settings["x_margin"], positions.max() + settings["x_margin"]),
           ylabel=y_label,
           xlabel=settings["x_label"])
    labelled_positions = set(positions[::settings["position_label_step"][variant_class]])
    ax.set_xticks(positions, [settings["position_label_format"][variant_class].format(
        wt=row.wt_aa, position=row.codon_position if row.codon_position in labelled_positions else "")
        for row in residues.itertuples()])
    ax.yaxis.set_major_locator(MultipleLocator(settings["y_tick_step"]))
    ax.yaxis.labelpad = settings["y_label_pad"]
    ax.yaxis.label.set_y(settings["y_label_y"])
    figure.canvas.draw()
    # Preserve panel geometry across exports and center outside legends on the axes.
    if not inset_legend:
        legend_right = 1 - figure.get_layout_engine().get()["w_pad"] / figure.get_figwidth()
        legend.set_bbox_to_anchor((legend_right, ax.get_position().y0 + ax.get_position().height / 2))
    figure.set_layout_engine("none")
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    labels = [ax.xaxis.label, ax.yaxis.label, legend, *ax.texts, *ax.get_xticklabels()]
    for star in [text for text in ax.texts if text.get_text() in dict(settings["star_thresholds"]).values()]:
        box = star.get_window_extent(renderer)
        require(box.y0 >= ax.get_window_extent(renderer).y0, "A significance star falls below its panel")
    labels += [label for value, label in zip(ax.get_yticks(), ax.get_yticklabels()) if 0 <= value <= maximum]
    boxes = [label.get_window_extent(renderer) for label in labels]
    require(all(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= figure.bbox.width and box.y1 <= figure.bbox.height
                for box in boxes), "A position profile label extends outside the page")
    require(all(not left.overlaps(right) for index, left in enumerate(boxes) for right in boxes[index + 1:]),
            "Position profile labels overlap")
    save_figure(figure, output)

    print(f"{output}, {sample_id}: {', '.join(sites)} plotted separately; "
          f"equal-weight arithmetic mean across {', '.join(classes)} variants at each position and site. "
          f"Linkage support >= {FIGURE_MIN_LINKAGE_SUPPORT}; zero-edit variants included.")
    print(profile.to_string(index=False, float_format="%.6f"))
    marked = profile.loc[profile.significant]
    print(f"Tested against the pooled synonymous wild type at a barcode-level dispersion of "
          + "; ".join(f"{site}={wild_type.loc[site].dispersion:.4f}" for site in sites)
          + f". Benjamini-Hochberg q < {SIGNIFICANCE_LEVEL:g} at "
          + "; ".join(f"{site}: {', '.join(map(str, sorted(marked.loc[marked.site_id.eq(site), 'codon_position'])))}"
                      for site in sites) + ".")
    print("Starred, both sites agreeing, at the least favourable of their two rates: "
          + ", ".join(f"{position} {symbol} (q={marks[position]:.4f})" for position, symbol in stars.items())
          + ". *** is q < 0.001, ** q < 0.01, and * q < 0.05. A position significant at one site alone "
            "carries no star: " + ", ".join(str(position) for position in
                                            sorted(set(marked.codon_position) - set(stars))) + ".")
    print("Pooled synonymous WT references (%): " + "; ".join(f"{site}={wt_percent[site]:.6f}" for site in sites))
    print("WT is not counted again at every position. No averaging across sites, pseudocount, "
          "WT normalization, or log-rate averaging.")
    print(f"Undefined site rates, omitted from position means: "
          f"{', '.join(cells.loc[cells.percent.isna(), 'variant'].drop_duplicates())}")
    print("The two sites share reads: their separation measures site-to-site variation, not replicate error. "
          "All derived summaries remain in memory and this log.")
