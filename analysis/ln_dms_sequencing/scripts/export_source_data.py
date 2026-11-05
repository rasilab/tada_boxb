"""Export the plotted values of each λN DMS manuscript panel as one tidy CSV.

Every number is taken from the same pooling helpers the plotting scripts call, so a panel and its
source data cannot disagree. Each percentage carries the edited and site-callable read counts it
was divided from, and each row carries the λN position, residues or codons, variant class,
recorder site, and barcode support that identify it. Uncovered measurements stay blank rather than
becoming zero. Run from the repository root.
"""

import os
import numpy as np
import pandas as pd

from analysis_helpers import (CLASS_MIN_BARCODES, FIGURE_MIN_LINKAGE_SUPPORT, MIN_CALLABLE_READS,
                              POOLED_WT_LABEL, annotation, comparison_sites, figure_sample, inserts,
                              measured_sites, pooled_variant_sites, require)

SOURCE_DATA_DIRECTORY = "source_data"
# The manuscript panel each figure is exported as.
PANELS = {
    "variant_read_coverage": "figure_6b",
    "recorder_site_profile": "figure_6c",
    "amino_acid_site_correlation_A2_A15": "figure_6d",
    "variant_class_distribution": "figure_6e",
    "position_editing_profile": "figure_6f",
    "amino_acid_editing_logo": "figure_6h",
}
PROFILE_CLASSES = {"position_editing_profile": "missense"}
DISTRIBUTION_CLASSES = ["synonymous", "missense", "stop"]
LOGO_CLASSES = ["missense"]
IDENTITY_NAMES = {"codon_position": "lambda_n_position", "wt_aa": "wild_type_residue",
                  "mut_aa": "substituted_residue", "wt_codon": "wild_type_codon",
                  "mut_codon": "substituted_codon"}
FLOAT_FORMAT = "%.6f"


def write_panel(frame, columns, figure_name):
    """Write one panel's plotted values under the manuscript panel name."""
    frame = frame.rename(columns=IDENTITY_NAMES)[columns]
    path = f"{SOURCE_DATA_DIRECTORY}/{PANELS[figure_name]}.csv"
    frame.to_csv(path, index=False, float_format=FLOAT_FORMAT, lineterminator="\n")
    print(f"{path} from {figure_name}: {len(frame):,} rows; columns {', '.join(frame.columns)}")


def site_values(pooled, site):
    """Barcode support, edited and callable reads, and the percentage, for one recorder site."""
    values = pooled.loc[pooled.site_id.eq(site),
                        ["variant", "barcodes", "edited", "callable", "fraction"]].copy()
    values[f"percent_editing_{site}"] = 100 * values.fraction
    return values.drop(columns="fraction").rename(columns={
        "barcodes": f"barcodes_{site}", "edited": f"edited_reads_{site}",
        "callable": f"callable_reads_{site}"})


sample_id = figure_sample()
comparison = comparison_sites()
designs = inserts()
os.makedirs(SOURCE_DATA_DIRECTORY, exist_ok=True)

recorder = measured_sites(primary_only=True).sort_values("rna_position")
require(recorder.feature.nunique() == 1, "The recorder profile requires one recorder feature")
require(set(comparison) <= set(recorder.site_id), "A comparison site is not a recorder site")
# One pooling pass per variant resolution; the two comparison sites are a slice of the recorder set.
pooled_recorder = pooled_variant_sites(sample_id, list(recorder.site_id), POOLED_WT_LABEL,
                                       FIGURE_MIN_LINKAGE_SUPPORT, variant_resolution="amino_acid")
pooled_codon = pooled_variant_sites(sample_id, comparison, POOLED_WT_LABEL,
                                    FIGURE_MIN_LINKAGE_SUPPORT, variant_resolution="codon")
pooled_amino = pooled_recorder.loc[pooled_recorder.site_id.isin(comparison)]

# Variant labels exactly as the pooling helper builds them, carrying their designed identities.
codon_identity = designs[["row_num", "codon_position", "wt_codon", "mut_codon", "wt_aa", "mut_aa",
                          "variant_class", "is_wild_type"]].copy()
codon_identity["variant"] = (codon_identity.wt_codon + codon_identity.codon_position.astype(str)
                             + codon_identity.mut_codon)
codon_identity.loc[codon_identity.is_wild_type, "variant"] = POOLED_WT_LABEL
amino_identity = designs[["codon_position", "wt_aa", "mut_aa", "variant_class"]].drop_duplicates().copy()
amino_identity["variant"] = (amino_identity.wt_aa + amino_identity.codon_position.astype(str)
                             + amino_identity.mut_aa)
amino_identity["codon_position"] = amino_identity.codon_position.astype("Int64")
# Synonymous designs pool into one wild-type group, which belongs to no single position or residue.
synonymous = amino_identity.variant_class.eq("synonymous")
for column, value in [("codon_position", pd.NA), ("wt_aa", pd.NA), ("mut_aa", pd.NA),
                      ("variant", POOLED_WT_LABEL)]:
    amino_identity.loc[synonymous, column] = value
amino_identity = amino_identity.drop_duplicates("variant").reset_index(drop=True)
require(set(amino_identity.variant) == set(pooled_amino.variant)
        and set(codon_identity.variant) == set(pooled_codon.variant),
        "Designed variant labels disagree with the pooled counts")

# Figure 6b: accepted reads per designed DNA insert, ranked. Reads do not depend on the site.
reads = pooled_codon.pivot(index="variant", columns="site_id", values="reads")[comparison]
require(reads.nunique(axis=1).eq(1).all(), "Accepted reads per variant disagree between sites")
coverage = (codon_identity.merge(reads[comparison[0]].rename("accepted_reads"), left_on="variant",
                                 right_index=True, validate="one_to_one")
            .sort_values(["accepted_reads", "variant"], ascending=[False, True]).reset_index(drop=True))
coverage.insert(0, "codon_variant_rank", np.arange(1, len(coverage) + 1))
write_panel(coverage.rename(columns={"variant": "codon_variant"}),
            ["codon_variant_rank", "codon_variant", "lambda_n_position", "wild_type_codon",
             "substituted_codon", "wild_type_residue", "substituted_residue", "variant_class",
             "accepted_reads"],
            "variant_read_coverage")

# Figure 6c: every amino-acid group at every recorder site, the values the boxes summarise.
profile = (pooled_recorder.merge(amino_identity.drop(columns="variant_class"), on="variant",
                                 validate="many_to_one")
           .merge(recorder[["site_id", "rna_position", "context"]], on="site_id", validate="many_to_one"))
profile["percent_editing"] = 100 * profile.fraction
profile = profile.sort_values(["rna_position", "codon_position", "variant"], na_position="first")
write_panel(profile.rename(columns={"site_id": "recorder_site", "rna_position": "recorder_rna_position",
                                    "context": "recorder_context", "variant": "amino_acid_variant",
                                    "edited": "edited_reads", "callable": "callable_reads"}),
            ["recorder_site", "recorder_rna_position", "recorder_context", "amino_acid_variant",
             "lambda_n_position", "wild_type_residue", "substituted_residue", "variant_class",
             "barcodes", "observed_codons", "edited_reads", "callable_reads", "percent_editing"],
            "recorder_site_profile")

# Figure 6d: one row per amino-acid group, the two coordinates of its scatter point.
correlation = amino_identity.copy()
for site in comparison:
    correlation = correlation.merge(site_values(pooled_amino, site), on="variant", validate="one_to_one")
correlation = correlation.sort_values(["codon_position", "variant"], na_position="first")
write_panel(correlation.rename(columns={"variant": "amino_acid_variant"}),
            ["amino_acid_variant", "lambda_n_position", "wild_type_residue", "substituted_residue",
             "variant_class"]
            + [f"{name}_{site}" for site in comparison
               for name in ["barcodes", "edited_reads", "callable_reads", "percent_editing"]],
            "amino_acid_site_correlation_A2_A15")

# Figure 6e: the codon variants that clear the barcode floor at both sites, and their site mean.
callable_codon = pooled_codon.pivot(index="variant", columns="site_id", values="callable")[comparison]
support_codon = pooled_codon.pivot(index="variant", columns="site_id", values="barcodes")[comparison]
rates_codon = pooled_codon.pivot(index="variant", columns="site_id", values="fraction")[comparison]
drawn = callable_codon.ge(MIN_CALLABLE_READS).all(axis=1) & support_codon.ge(CLASS_MIN_BARCODES).all(axis=1)
distribution = codon_identity.loc[codon_identity.variant.isin(drawn.index[drawn])].copy()
for site in comparison:
    distribution = distribution.merge(site_values(pooled_codon, site), on="variant", validate="one_to_one")
distribution["mean_percent_editing"] = (
    100 * rates_codon.loc[distribution.variant].mean(axis=1, skipna=False).to_numpy())
require(distribution.variant_class.isin(DISTRIBUTION_CLASSES).all()
        and distribution.mean_percent_editing.between(0, 100).all(),
        "A class-distribution variant has an unknown class or an impossible rate")
distribution = distribution.sort_values(["variant_class", "codon_position", "variant"])
write_panel(distribution.rename(columns={"variant": "codon_variant"}),
            ["codon_variant", "lambda_n_position", "wild_type_codon", "substituted_codon",
             "wild_type_residue", "substituted_residue", "variant_class"]
            + [f"{name}_{site}" for site in comparison
               for name in ["barcodes", "edited_reads", "callable_reads", "percent_editing"]]
            + ["mean_percent_editing"],
            "variant_class_distribution")

# Figure 6f: the equal-weight position means each profile draws, with their WT reference.
residues = designs[["codon_position", "wt_aa"]].drop_duplicates().sort_values("codon_position")
callable_amino = pooled_amino.pivot(index="variant", columns="site_id", values="callable")[comparison]
covered_amino = callable_amino.ge(MIN_CALLABLE_READS).all(axis=1)
reference = pooled_amino.loc[pooled_amino.variant.eq(POOLED_WT_LABEL)].set_index("site_id")
for figure_name, variant_class in PROFILE_CLASSES.items():
    cells = designs.loc[designs.variant_class.eq(variant_class),
                        ["codon_position", "wt_aa", "mut_aa"]].drop_duplicates().copy()
    cells["variant"] = cells.wt_aa + cells.codon_position.astype(str) + cells.mut_aa
    cells = cells.merge(pooled_amino[["variant", "site_id", "edited", "callable", "barcodes", "fraction"]],
                        on="variant", validate="one_to_many")
    # A group uncovered at either site leaves both site means, exactly as the panel draws them.
    cells["percent"] = 100 * cells.fraction.where(cells.variant.map(covered_amino))
    counted = cells.loc[cells.percent.notna()]
    summary = counted.groupby(["codon_position", "site_id"], as_index=False).agg(
        mean_percent_editing=("percent", "mean"), covered_variants=("percent", "size"),
        edited_reads=("edited", "sum"), callable_reads=("callable", "sum"), barcodes=("barcodes", "sum"))
    designed = (cells.groupby(["codon_position", "site_id"], as_index=False).size()
                .rename(columns={"size": "designed_variants"}))
    summary = (residues.merge(summary, on="codon_position", validate="one_to_many")
               .merge(designed, on=["codon_position", "site_id"], validate="one_to_one"))
    summary["variant_class"] = variant_class
    summary["wild_type_percent_editing"] = 100 * summary.site_id.map(reference.fraction)
    summary["wild_type_edited_reads"] = summary.site_id.map(reference.edited)
    summary["wild_type_callable_reads"] = summary.site_id.map(reference["callable"])
    require(len(summary) == len(residues) * len(comparison)
            and summary.mean_percent_editing.between(0, 100).all(),
            f"A {variant_class} position lacks a valid mean editing rate")
    summary["site_id"] = pd.Categorical(summary.site_id, categories=comparison, ordered=True)
    summary = summary.sort_values(["codon_position", "site_id"])
    write_panel(summary.rename(columns={"site_id": "recorder_site"}),
                ["lambda_n_position", "wild_type_residue", "recorder_site", "variant_class",
                 "mean_percent_editing", "covered_variants", "designed_variants", "barcodes",
                 "edited_reads", "callable_reads", "wild_type_percent_editing",
                 "wild_type_edited_reads", "wild_type_callable_reads"],
                figure_name)

# Figure 6h: the normalized letter heights of every stack, with the rates they are normalized from.
groups = annotation("amino_acid_groups")
logo_residues = designs.loc[designs.variant_class.isin([*LOGO_CLASSES, "synonymous"]), "mut_aa"].unique()
grid = residues.merge(groups.loc[groups.amino_acid.isin(logo_residues)], how="cross")
grid["is_wild_type_reference"] = grid.amino_acid.eq(grid.wt_aa)
grid["variant"] = grid.wt_aa + grid.codon_position.astype(str) + grid.amino_acid
grid.loc[grid.is_wild_type_reference, "variant"] = POOLED_WT_LABEL
for site in comparison:
    grid = grid.merge(site_values(pooled_amino, site), on="variant", validate="many_to_one")
grid["mean_percent_editing"] = grid[[f"percent_editing_{site}" for site in comparison]].mean(
    axis=1, skipna=False)
grid["relative_letter_height"] = (grid.mean_percent_editing
                                  / grid.codon_position.map(grid.groupby("codon_position")
                                                            .mean_percent_editing.sum(min_count=1)))
require(np.allclose(grid.groupby("codon_position").relative_letter_height.sum(), 1),
        "A logo stack does not reach unit height")
grid = grid.sort_values(["codon_position", "relative_letter_height"])
write_panel(grid.rename(columns={"amino_acid": "logo_residue", "group": "amino_acid_group"}),
            ["lambda_n_position", "wild_type_residue", "logo_residue", "amino_acid_group",
             "is_wild_type_reference"]
            + [f"{name}_{site}" for site in comparison
               for name in ["barcodes", "edited_reads", "callable_reads", "percent_editing"]]
            + ["mean_percent_editing", "relative_letter_height"],
            "amino_acid_editing_logo")

print(f"{sample_id}: {len(PANELS)} panel tables written to {SOURCE_DATA_DIRECTORY}/ at a linkage "
      f"floor of {FIGURE_MIN_LINKAGE_SUPPORT} designed reads, from the same pooled counts the "
      "figures draw. Blank percentages are uncovered measurements, not zero editing.")
