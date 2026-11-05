"""A8/A2 agreement using direct editing percentages and the in vitro scatter layout.

Each comparison names the sequences it draws, the 57 single substitutions or every mutant the
combinatorial library designs. Edited and callable counts are pooled across matching constructs
separately at each site. Linear axes include every drawn mutant, zero-edit mutants included; WT is not
plotted. Sites share molecules and are not independent replicates. An edit inside boxB moves a
molecule into a genotype further from wild type, which reaches the multi-mutants the designed set
adds more than the single substitutions.
"""

import numpy as np
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import gaussian_kde, pearsonr, spearmanr

from analysis_helpers import (annotation, boxb_wildtype, comparison_samples, designed_boxb_sequences,
                              excluded_constructs, require, sequence_effects,
                              substitution_measurements)
from plot_style import centre_shared_xlabel, palettes, place_shared_ylabel, plt, save_figure

FIGURE_SIZE = (1.50, 1.24)
# The designed set draws an order of magnitude more points into the same square, so its markers are
# smaller and translucent, which keeps the dense centre of the distribution readable. The in vitro
# all-variant panels take the same two styles.
VARIANT_POINT_STYLES = {"single": dict(s=7, alpha=1), "designed": dict(s=1.6, alpha=.45)}
TICK_BINS = 3
TICK_INTERVALS = [1, 2, 2.5, 5, 10]
REFERENCE_COLOR = "#999999"
REFERENCE_WIDTH = .4
CORRELATION_POSITION = (.04, .96)
DENSITY_POINTS = 256
DENSITY_WIDTH_PERCENT = 16
DENSITY_GAP_IN = .04
DENSITY_ALPHA = .6
# Match the in vitro plot's nrd0 bandwidth.
BANDWIDTH_FACTOR = .9
IQR_TO_SD = 1.34

labels, _, colors = palettes()
comparisons = annotation("boxb_substitution_correlations")
require(not comparisons.source.duplicated().any() and comparisons.notna().all().all(),
        "Incomplete or repeated single-mutant comparison")
require(comparisons.variant_set.isin(VARIANT_POINT_STYLES).all(), "Every comparison must name a known variant set")
require(excluded_constructs().empty, "Flagged placements require an explicit with/without comparison")
designed = designed_boxb_sequences()
mutants = designed.loc[designed.n_mutations.gt(0)].copy()
variant_sets = {"single": mutants.loc[mutants.n_mutations.eq(1)], "designed": mutants}
_, _, wildtype = boxb_wildtype()
mutation_names = {}
for sequence in mutants.hairpin:
    changes = [f"{reference}{index + 1}{base}".replace("T", "U")
               for index, (reference, base) in enumerate(zip(wildtype, sequence)) if reference != base]
    require(len(changes) > 0, "A selected sequence carries no substitution")
    mutation_names[sequence] = "+".join(changes)

for comparison in comparisons.itertuples():
    variants = variant_sets[comparison.variant_set]
    selected = comparison_samples(comparison.comparison_set)
    selected = selected.loc[selected.sample_id.eq(comparison.sample_id)]
    require(len(selected) == 1, "A scatter must select one annotated library")
    sites = [comparison.x_position_id, comparison.y_position_id]
    measurements = substitution_measurements(comparison.comparison_set).set_index("position_id")
    require(len(set(sites)) == 2 and set(sites) <= set(measurements.index), "Invalid site comparison")
    measurements = measurements.loc[sites].reset_index()
    effects, constructs = sequence_effects(selected, measurements)
    effects = effects.merge(variants[["hairpin"]], on="hairpin", validate="many_to_one")
    require(len(effects) == len(selected) * len(sites) * len(variants)
            and not effects.duplicated(["sample_id", "position_id", "hairpin"]).any(),
            "A library is missing a designed mutant measurement")
    require(effects.pooled_umis.gt(0).all() and effects.pooled_fraction_edited.between(0, 1).all(),
            "A sequence lacks callable coverage or a valid editing fraction")
    retained = constructs.loc[constructs.hairpin.isin(variants.hairpin)]

    panels = []
    for library in selected.itertuples():
        rows = effects.loc[effects.sample_id.eq(library.sample_id)]
        points = 100 * rows.pivot(index="hairpin", columns="position_id", values="pooled_fraction_edited")
        drawn = points.reindex(variants.hairpin)
        require(np.isfinite(drawn[sites]).all().all(), "A mutant editing percentage is missing")
        x, y = drawn[sites[0]].to_numpy(), drawn[sites[1]].to_numpy()
        require(np.std(x) > 0 and np.std(y) > 0, "A constant measurement has no correlation")
        correlation = pearsonr(x, y)
        panels.append((library, x, y, correlation.statistic))
        print(f"{comparison.source} {library.sample_id} {labels[library.enzyme]}: n={len(x)}/{len(variants)}, "
              f"Pearson editing-percent r={correlation.statistic:.6f}, p={correlation.pvalue:.3g}; "
              f"Spearman rho={spearmanr(x, y).statistic:.6f}; "
              f"SD in percentage points {sites[0]}={np.std(x, ddof=1):.6f}, {sites[1]}={np.std(y, ddof=1):.6f}")
        zeros = rows.loc[rows.hairpin.isin(variants.hairpin) & rows.pooled_edited_umis.eq(0)].copy()
        if len(zeros):
            zeros["mutation"] = zeros.hairpin.map(mutation_names)
            print("  Included on linear axes at zero editing:")
            print(zeros[["mutation", "position_id", "pooled_edited_umis", "pooled_umis"]]
                  .sort_values(["mutation", "position_id"]).to_string(index=False))

    values = np.concatenate([array for _, x, y, _ in panels for array in [x, y]])
    ticks = MaxNLocator(nbins=TICK_BINS, steps=TICK_INTERVALS).tick_values(0, values.max())
    limits = (0, ticks[-1])
    tick_labels = [f"{value:.0f}" for value in ticks]
    figure, axes = plt.subplots(1, len(panels), figsize=FIGURE_SIZE, layout="constrained",
                                sharex=True, sharey=True, squeeze=False)
    for ax, (library, x, y, correlation) in zip(axes.flat, panels):
        density_ax = make_axes_locatable(ax).append_axes(
            "right", size=f"{DENSITY_WIDTH_PERCENT}%", pad=DENSITY_GAP_IN, sharey=ax)
        ax.plot(limits, limits, linestyle="--", color=REFERENCE_COLOR, linewidth=REFERENCE_WIDTH, zorder=0)
        ax.scatter(x, y, color=colors[library.enzyme], linewidths=0, clip_on=False,
                   **VARIANT_POINT_STYLES[comparison.variant_set])
        bandwidth = BANDWIDTH_FACTOR * min(np.std(y, ddof=1), np.subtract(*np.percentile(y, [75, 25])) / IQR_TO_SD) * len(y) ** -.2
        require(bandwidth > 0, "The marginal distribution has zero bandwidth")
        density_y = np.linspace(*limits, DENSITY_POINTS)
        density = gaussian_kde(y, bw_method=bandwidth / np.std(y, ddof=1))(density_y)
        density_ax.fill_betweenx(density_y, 0, density / density.max(), color=colors[library.enzyme],
                                alpha=DENSITY_ALPHA, linewidth=0)
        density_ax.axvline(0, color=REFERENCE_COLOR, linewidth=REFERENCE_WIDTH)
        density_ax.set_xlim(0, 1)
        density_ax.axis("off")
        count_label = f"{len(x)}/{len(variants)}" if len(x) != len(variants) else str(len(x))
        ax.text(*CORRELATION_POSITION, f"r = {correlation:.2f}\nn = {count_label}",
                transform=ax.transAxes, va="top")
        ax.set(xlim=limits, ylim=limits,
               xticks=ticks, yticks=ticks, xticklabels=tick_labels, yticklabels=tick_labels)
        ax.set_aspect("equal")
        if ax is not axes[0, 0]:
            ax.tick_params(left=False, labelleft=False)
            ax.spines["left"].set_visible(False)
    shared_x = figure.supxlabel(f"Editing at {sites[0]} (%)")
    place_shared_ylabel(figure, axes[0, 0], f"Editing at {sites[1]} (%)")
    centre_shared_xlabel(figure, shared_x)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    for left, right in zip(axes.flat, list(axes.flat)[1:]):
        require(max(text.get_window_extent(renderer).x1 for text in left.get_xticklabels())
                <= min(text.get_window_extent(renderer).x0 for text in right.get_xticklabels()),
                "Neighbouring panels' tick labels overlap")
    save_figure(figure, comparison.source)
    print("No in vivo placements are excluded; with/without-exclusion estimates are identical.")
