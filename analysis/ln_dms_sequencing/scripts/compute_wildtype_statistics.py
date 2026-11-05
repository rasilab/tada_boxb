"""Derive the wild-type and robustness statistics this analysis reports outside a figure.

Run from the repository root. Three groups of numbers are computed from the same linkage and count
tables the figures use, through the same helpers: how wild-type barcode editing depends on read
depth, the calibration of the barcode test against the read-level test, and the read-floor and
barcode-floor robustness series. One tidy row is written per statistic.

A pseudo-variant is tested against a reference rebuilt from the wild-type barcodes it did not draw,
because that matches how the test is applied to a real variant, whose barcodes are never wild-type
barcodes, and a calibration should simulate the procedure being calibrated. Leaving the drawn
barcodes in the reference instead correlates the two rates and makes the barcode test mildly
conservative. Measured over twelve seeds, holding out gives barcode rates of 5.1% at A2 and 5.8% at
A15 against 4.1% and 5.0% with the reference left whole, and leaves the Fisher rates alone, 88.8%
and 89.4% against 88.4% and 89.2%. The calibration draws are random and use the seed recorded below.
"""

import os

import numpy as np
import pandas as pd

from analysis_helpers import (FIGURE_MIN_LINKAGE_SUPPORT, MIN_CALLABLE_READS, POOLED_WT_LABEL,
                              SIGNIFICANCE_LEVEL, barcode_dispersion, clustered_rate_error,
                              clustered_rate_test, comparison_sites, figure_sample, inserts,
                              pooled_variant_sites, read_rate_test, require, variant_barcode_sites,
                              wild_type_comparisons, wild_type_reference)

OUTPUT_TABLE = "analysis/ln_dms_sequencing/tables/wildtype_statistics.csv"
# Wild-type barcodes thinner than this carry too few reads to say anything about their own rate,
# so the dispersion and the depth contrast are both read off the barcodes above it.
DISPERSION_MIN_CALLABLE = 50
# The two depth bins contrasted here, in accepted RNA reads per barcode, and the editing
# percentage that counts as "edits below 1%".
DEPTH_LOW_READS = 200
DEPTH_HIGH_READS = 1000
LOW_EDITING_PERCENT = 1.0
# Calibration: pseudo-variants drawn from the wild-type barcodes alone, sized like the real
# amino-acid variants, and tested at the level the figures use.
NOMINAL_LEVEL = SIGNIFICANCE_LEVEL
CALIBRATION_DRAWS = 400
CALIBRATION_SEED = 20261029
# Robustness series: callable reads required of a variant, accepted reads required of a barcode,
# and barcodes required of a codon variant at both sites.
VARIANT_READ_FLOORS = [1, 100, 200, 500]
BARCODE_READ_FLOORS = [1, 50, 100, 200, 500]
CODON_BARCODE_FLOORS = [1, 3, 5]
ARGININE_RESIDUE = "R"
PERCENT_QUANTILE = 5

sample_id = figure_sample()
sites = comparison_sites()
designs = inserts()
residues = designs[["codon_position", "wt_aa"]].drop_duplicates().sort_values("codon_position")
arginine_positions = sorted(residues.loc[residues.wt_aa.eq(ARGININE_RESIDUE), "codon_position"])
require(len(arginine_positions) > 0, "The λN segment carries no arginine position")

barcode_sites, _, _ = variant_barcode_sites(sample_id, sites, POOLED_WT_LABEL,
                                            FIGURE_MIN_LINKAGE_SUPPORT)
references = wild_type_reference(barcode_sites, POOLED_WT_LABEL).set_index("site_id")
wild_type = barcode_sites.loc[barcode_sites.variant.eq(POOLED_WT_LABEL)]
rows = []


def record(statistic, value, unit, site="both", threshold="", n=np.nan, detail=""):
    rows.append(dict(statistic=statistic, site=site, threshold=threshold, value=value, unit=unit,
                     n=n, detail=detail))


# Wild-type barcodes: how far they disagree, and how that disagreement tracks read depth.
for site in sites:
    reference = references.loc[site]
    record("wild_type_pooled_rate", 100 * reference.reference_rate, "percent", site=site,
           n=reference.reference_barcodes,
           detail=f"{int(reference.reference_edited)} edited of "
                  f"{int(reference.reference_callable)} callable")
    record("wild_type_intraclass_correlation", reference.dispersion, "correlation", site=site,
           n=reference.reference_barcodes, detail="every covered wild-type barcode")

deep = wild_type.groupby("barcode")["callable"].min().ge(DISPERSION_MIN_CALLABLE)
deep_barcodes = set(deep.loc[deep].index)
require(len(deep_barcodes) > 1, "Too few wild-type barcodes carry the dispersion floor")
for site in sites:
    reference = references.loc[site]
    frame = wild_type.loc[wild_type.site_id.eq(site) & wild_type.barcode.isin(deep_barcodes)]
    expected = frame["callable"] * reference.reference_rate
    pearson = (((frame.edited - expected) ** 2)
               / (expected * (1 - reference.reference_rate))).sum() / (len(frame) - 1)
    record("wild_type_pearson_dispersion", pearson, "ratio", site=site,
           threshold=f">= {DISPERSION_MIN_CALLABLE} callable reads", n=len(frame),
           detail="at the pooled wild-type rate of this site")
    record("binomial_error_understatement", np.sqrt(pearson), "fold", site=site,
           threshold=f">= {DISPERSION_MIN_CALLABLE} callable reads", n=len(frame),
           detail="square root of the Pearson dispersion")
    record("wild_type_median_barcode_rate", 100 * (frame.edited / frame["callable"]).median(),
           "percent", site=site, threshold=f">= {DISPERSION_MIN_CALLABLE} callable reads",
           n=len(frame), detail="median of the per-barcode rates")

# The depth contrast pools the two sites within each barcode, because a barcode is one clone
# measured twice, and bins it on the accepted RNA reads that clone carries.
depth = (wild_type.loc[wild_type.barcode.isin(deep_barcodes)].groupby("barcode")
         .agg(edited=("edited", "sum"), callable=("callable", "sum"), reads=("reads", "max")))
depth["below_threshold"] = 100 * depth.edited / depth["callable"] < LOW_EDITING_PERCENT
for label, selected, bound in [("below", depth.reads.lt(DEPTH_LOW_READS), DEPTH_LOW_READS),
                               ("above", depth.reads.gt(DEPTH_HIGH_READS), DEPTH_HIGH_READS)]:
    group = depth.loc[selected]
    require(len(group) > 0, f"No wild-type barcode falls {label} {bound} reads")
    record(f"wild_type_share_below_{LOW_EDITING_PERCENT:g}_percent_{label}_{bound}_reads",
           100 * group.below_threshold.mean(), "percent",
           threshold=f">= {DISPERSION_MIN_CALLABLE} callable reads at both sites", n=len(group),
           detail=f"{int(group.below_threshold.sum())} of {len(group)} barcodes, sites summed")

# Calibration: wild type against itself, at the sizes the real variants carry.
generator = np.random.default_rng(CALIBRATION_SEED)
for site in sites:
    covered = barcode_sites.loc[barcode_sites.site_id.eq(site) & barcode_sites["callable"].gt(0)]
    pool = covered.loc[covered.variant.eq(POOLED_WT_LABEL)]
    sizes = covered.loc[covered.variant.ne(POOLED_WT_LABEL)].groupby("variant").size().to_numpy()
    require(len(sizes) > 0 and len(pool) > 1, "Calibration needs real variants and a wild-type pool")
    require(sizes.max() < len(pool) - 1, "A pseudo-variant would empty the wild-type reference")
    edited = pool.edited.to_numpy(float)
    callable_counts = pool["callable"].to_numpy(float)
    rejected = {"barcode": 0, "read": 0}
    for _ in range(CALIBRATION_DRAWS):
        size = int(generator.choice(sizes))
        drawn = generator.choice(len(pool), size=size, replace=False)
        # A real variant's barcodes are never wild-type barcodes, so the reference a pseudo-variant
        # is tested against is rebuilt from the barcodes it did not take. Leaving them in would
        # correlate the two rates and calibrate a comparison the figures never make.
        rest = np.setdiff1d(np.arange(len(pool)), drawn)
        rest_rate = edited[rest].sum() / callable_counts[rest].sum()
        rest_dispersion = barcode_dispersion(edited[rest], callable_counts[rest], rest_rate)
        rest_error = clustered_rate_error(edited[rest], callable_counts[rest], rest_rate,
                                          rest_dispersion)
        _, _, p_barcode = clustered_rate_test(edited[drawn], callable_counts[drawn], rest_rate,
                                              rest_error, rest_dispersion)
        p_read = read_rate_test(int(edited[drawn].sum()), int(callable_counts[drawn].sum()),
                                int(edited[rest].sum()), int(callable_counts[rest].sum()))
        rejected["barcode"] += p_barcode < NOMINAL_LEVEL
        rejected["read"] += p_read < NOMINAL_LEVEL
    for test, count in rejected.items():
        record(f"calibration_{test}_test_rejection", 100 * count / CALIBRATION_DRAWS, "percent",
               site=site, threshold=f"nominal {NOMINAL_LEVEL:g}", n=CALIBRATION_DRAWS,
               detail=f"{count} of {CALIBRATION_DRAWS} draws, seed {CALIBRATION_SEED}, "
                      f"pseudo-variant sizes {sizes.min()}-{sizes.max()} barcodes, "
                      f"reference held out")

comparisons, _ = wild_type_comparisons(sample_id, sites, POOLED_WT_LABEL, FIGURE_MIN_LINKAGE_SUPPORT)
tested = comparisons.loc[comparisons.variant.ne(POOLED_WT_LABEL)]
barcode_significant = tested.q_value_barcode.lt(SIGNIFICANCE_LEVEL)
read_significant = tested.q_value_read.lt(SIGNIFICANCE_LEVEL)
both_sites = (tested.loc[barcode_significant].groupby("variant").size().ge(len(sites)).sum())
record("variant_site_comparisons", len(tested), "comparisons", n=len(tested),
       detail=f"{len(sites)} sites x the amino-acid variants covered at each")
record("barcode_test_significant", int(barcode_significant.sum()), "comparisons", n=len(tested),
       threshold=f"q < {SIGNIFICANCE_LEVEL:g}", detail="the test the figures mark")
record("read_test_significant", int(read_significant.sum()), "comparisons", n=len(tested),
       threshold=f"q < {SIGNIFICANCE_LEVEL:g}", detail="Fisher on pooled reads, off the figures")
record("barcode_test_significant_both_sites", int(both_sites), "variants", n=len(tested),
       threshold=f"q < {SIGNIFICANCE_LEVEL:g}", detail="variants reaching it at both sites")

# Read depth does not make the low missense tail: neither a floor on a variant's reads nor one on
# a barcode's reads moves the missense rates.
pooled = pooled_variant_sites(sample_id, sites, POOLED_WT_LABEL, FIGURE_MIN_LINKAGE_SUPPORT)
missense = pooled.loc[pooled.variant_class.eq("missense")]
coverage = missense.pivot(index="variant", columns="site_id", values="callable")[sites]
rates = missense.pivot(index="variant", columns="site_id", values="fraction")[sites]
covered = coverage.ge(MIN_CALLABLE_READS).all(axis=1)
require(covered.sum() > 0, "No missense amino-acid group is covered at both sites")
for site in sites:
    callable_reads = coverage.loc[covered, site]
    zero = rates.loc[covered, site].eq(0)
    record("missense_median_callable_reads", callable_reads.median(), "reads", site=site,
           n=int(covered.sum()), detail="covered missense amino-acid groups")
    record("missense_fifth_percentile_callable_reads",
           np.percentile(callable_reads, PERCENT_QUANTILE), "reads", site=site,
           n=int(covered.sum()), detail=f"{PERCENT_QUANTILE}th percentile over the same groups")
    record("missense_zero_groups", int(zero.sum()), "groups", site=site, n=int(covered.sum()),
           detail="groups reading exactly 0% at this site")
    record("missense_zero_group_median_callable_reads", callable_reads[zero].median(), "reads",
           site=site, n=int(zero.sum()), detail="callable reads those zero groups carry")
    record("missense_expected_zero_groups",
           ((1 - references.loc[site].reference_rate) ** callable_reads).sum(), "groups", site=site,
           n=int(covered.sum()), detail="expected under the pooled wild-type rate")

for site in sites:
    means = []
    for floor in VARIANT_READ_FLOORS:
        kept = coverage.ge(floor).all(axis=1)
        mean_percent = 100 * rates.loc[kept, site].mean()
        means.append(mean_percent)
        record("missense_mean_rate_at_variant_read_floor", mean_percent, "percent", site=site,
               threshold=f">= {floor} callable reads", n=int(kept.sum()),
               detail="equal weight over the covered missense groups retained")
    record("missense_mean_rate_largest_shift", max(means) - min(means), "percentage points",
           site=site, threshold=f"floors {', '.join(map(str, VARIANT_READ_FLOORS))}",
           n=len(VARIANT_READ_FLOORS), detail="largest gap across the read floors above")

cells = designs.loc[designs.variant_class.eq("missense"),
                    ["codon_position", "wt_aa", "mut_aa"]].drop_duplicates().copy()
cells["variant"] = cells.wt_aa + cells.codon_position.astype(str) + cells.mut_aa
for floor in BARCODE_READ_FLOORS:
    kept = barcode_sites.loc[barcode_sites.reads.ge(floor)]
    summary = kept.groupby(["variant", "site_id"]).agg(
        edited=("edited", "sum"), callable=("callable", "sum")).reset_index()
    summary["percent"] = 100 * summary.edited / summary["callable"].replace(0, np.nan)
    frame = cells.merge(summary, on="variant", how="left")
    for site in sites:
        position = frame.loc[frame.site_id.eq(site)].groupby("codon_position").percent.mean()
        is_arginine = position.index.isin(arginine_positions)
        record("missense_arginine_position_mean", position[is_arginine].mean(), "percent",
               site=site, threshold=f">= {floor} reads per barcode", n=int(is_arginine.sum()),
               detail="equal weight over the arginine positions")
        record("missense_other_position_mean", position[~is_arginine].mean(), "percent", site=site,
               threshold=f">= {floor} reads per barcode", n=int((~is_arginine).sum()),
               detail="equal weight over the remaining positions")

# What the low tail does track is barcode support, at the codon resolution the class panel draws.
codons = pooled_variant_sites(sample_id, sites, POOLED_WT_LABEL, FIGURE_MIN_LINKAGE_SUPPORT,
                              variant_resolution="codon")
codon_coverage = codons.pivot(index="variant", columns="site_id", values="callable")[sites]
codon_rates = codons.pivot(index="variant", columns="site_id", values="fraction")[sites]
codon_support = codons.pivot(index="variant", columns="site_id", values="barcodes")[sites]
codon_class = codons[["variant", "variant_class"]].drop_duplicates().set_index("variant").variant_class
codon_percent = 100 * codon_rates.mean(axis=1, skipna=False)
for floor in CODON_BARCODE_FLOORS:
    retained = codon_coverage.ge(MIN_CALLABLE_READS).all(axis=1) & codon_support.ge(floor).all(axis=1)
    for variant_class in ["missense", "synonymous"]:
        # Zero means have no spread to report on a ratio scale, so the series follows the positive
        # variants the class panel draws, and counts the zeros it leaves off.
        selected = codon_percent.loc[retained & codon_class.eq(variant_class)]
        positive = selected.loc[selected.gt(0)]
        require(len(positive) > 1, f"Too few positive {variant_class} codon variants at floor {floor}")
        first, third = np.percentile(positive, [25, 75])
        record(f"{variant_class}_codon_variants_retained", len(positive), "variants",
               threshold=f">= {floor} barcodes", n=len(selected),
               detail=f"positive mean editing, of {len(selected)} covered at both sites")
        record(f"{variant_class}_codon_interquartile_spread", third / first, "fold",
               threshold=f">= {floor} barcodes", n=len(positive),
               detail="75th over 25th percentile of mean A2/A15 editing")
        record(f"{variant_class}_codon_share_below_{LOW_EDITING_PERCENT:g}_percent",
               100 * positive.lt(LOW_EDITING_PERCENT).mean(), "percent",
               threshold=f">= {floor} barcodes", n=len(positive),
               detail=f"{int(positive.lt(LOW_EDITING_PERCENT).sum())} of {len(positive)} variants")
        record(f"{variant_class}_codon_median_rate", positive.median(), "percent",
               threshold=f">= {floor} barcodes", n=len(positive), detail="median mean A2/A15 editing")

table = pd.DataFrame(rows)
require(table.value.notna().all(), "A statistic is undefined")
os.makedirs(os.path.dirname(OUTPUT_TABLE), exist_ok=True)
table.to_csv(OUTPUT_TABLE, index=False, float_format="%.6g")
print(f"{sample_id}: {len(table)} wild-type and robustness statistics over {', '.join(sites)}, "
      f"written to {OUTPUT_TABLE}.")
print(f"Linkage support >= {FIGURE_MIN_LINKAGE_SUPPORT}; wild-type barcodes "
      f"{int(references.reference_barcodes.iloc[0])}, of which {len(deep_barcodes)} carry "
      f"{DISPERSION_MIN_CALLABLE} callable reads at both sites.")
print(table.to_string(index=False, float_format="%.6g"))
