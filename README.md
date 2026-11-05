# Deaminase-based RNA recording enables high throughput mutational profiling of protein-RNA interactions

Rachael A. Bakker, Oliver B. Nicholson, Heungwon Park, Yu-Lan Xiao, Weixin Tang,
[Arvind Rasi Subramaniam](https://github.com/rasi), and
[Christopher P. Lapointe](https://github.com/cplapointe)

| | |
| --- | --- |
| Preprint | bioRxiv [10.1101/2025.04.11.648485](https://www.biorxiv.org/content/10.1101/2025.04.11.648485) |
| Sequencing data | NCBI BioProject [PRJNA1255650](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1255650) |
| Per-run accessions | [sra_accessions.csv](./data/sra_accessions.csv) |
| Plasmid sequences | [data/plasmid_maps](./data/plasmid_maps) |
| Per-panel data | [source_data](./source_data) |
| Questions and reagent requests | [open an issue](https://github.com/rasilab/tada_boxb/issues/new/choose) |
| Licence | [MIT](./LICENSE) for the code, [CC BY 4.0](./LICENSE-DATA) for the data |

## Figures

| Panel | Figure | Code | Source data |
| --- | --- | --- | --- |
| 1C | [reporter_editing](./analysis/nanopore_reporter_editing/figures/reporter_editing.pdf) | [plot_reporter_editing.py](./analysis/nanopore_reporter_editing/scripts/plot_reporter_editing.py) | [figure_1c.csv](./source_data/figure_1c.csv) |
| 2C | [recorder_edit_counts_enzyme_concentration](./analysis/boxb_in_vitro_sequencing/figures/recorder_edit_counts_enzyme_concentration.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_2c.csv](./source_data/figure_2c.csv) |
| 2D | [recorder_site_distance_tada_lambdaN_250nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_distance_tada_lambdaN_250nM.pdf) <br> [recorder_site_distance_tada_250nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_distance_tada_250nM.pdf) | [plot_recorder_site_distances.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_site_distances.py) | [figure_2d.csv](./source_data/figure_2d.csv) |
| 2E | [recorder_context_editing_A8](./analysis/boxb_in_vitro_sequencing/figures/recorder_context_editing_A8.pdf) | [plot_recorder_contexts.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_contexts.py) | [figure_2e.csv](./source_data/figure_2e.csv) |
| 3B | [boxb_all_variant_A8_vs_A2_250nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_all_variant_A8_vs_A2_250nM.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_3b.csv](./source_data/figure_3b.csv) |
| 3C | [boxb_mutation_distance_editing_mean_A8_A2_250nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_mutation_distance_editing_mean_A8_A2_250nM.pdf) | [plot_boxb_mutation_distance.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_mutation_distance.py) | [figure_3c.csv](./source_data/figure_3c.csv) |
| 3D | [boxb_single_mutant_editing_mean_A8_A2_250nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_single_mutant_editing_mean_A8_A2_250nM.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_3d.csv](./source_data/figure_3d.csv) |
| 3E | [boxb_lambdan_structure](./analysis/boxb_structure_views/figures/boxb_lambdan_structure.pdf) | rendered from the Legault et al. 1998 coordinates, obtained from that paper's corresponding author | — |
| 3F | [boxb_multimutant_heatmaps_mean_A8_A2_lambdaN_250nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_multimutant_heatmaps_mean_A8_A2_lambdaN_250nM.pdf) | [plot_boxb_multimutant_heatmaps.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_multimutant_heatmaps.py) | [figure_3f.csv](./source_data/figure_3f.csv) |
| 4B | [recorder_edit_counts_pag_vs_nanobody_250nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_edit_counts_pag_vs_nanobody_250nM.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_4b.csv](./source_data/figure_4b.csv) |
| 4C | [recorder_site_editing_pag_components_A8_250nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_editing_pag_components_A8_250nM.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_4c.csv](./source_data/figure_4c.csv) |
| 4D | [recorder_site_editing_pag_conditions_A8](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_editing_pag_conditions_A8.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_4d.csv](./source_data/figure_4d.csv) |
| 4E | [boxb_all_variant_A8_vs_A2_editors](./analysis/boxb_in_vitro_sequencing/figures/boxb_all_variant_A8_vs_A2_editors.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_4e.csv](./source_data/figure_4e.csv) |
| 4F | [boxb_mutation_distance_editing_mean_A8_A2_editors](./analysis/boxb_in_vitro_sequencing/figures/boxb_mutation_distance_editing_mean_A8_A2_editors.pdf) | [plot_boxb_mutation_distance.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_mutation_distance.py) | [figure_4f.csv](./source_data/figure_4f.csv) |
| 4G | [boxb_single_mutant_editing_mean_A8_A2_editors](./analysis/boxb_in_vitro_sequencing/figures/boxb_single_mutant_editing_mean_A8_A2_editors.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_4g.csv](./source_data/figure_4g.csv) |
| 5B | [boxb_read_coverage_tada_lambdaN_original](./analysis/boxb_in_vivo_sequencing/figures/boxb_read_coverage_tada_lambdaN_original.pdf) | [plot_boxb_read_coverage.py](./analysis/boxb_in_vivo_sequencing/scripts/plot_boxb_read_coverage.py) | [figure_5b.csv](./source_data/figure_5b.csv) |
| 5C | [recorder_edit_counts_editors_original](./analysis/boxb_in_vivo_sequencing/figures/recorder_edit_counts_editors_original.pdf) | [plot_recorder_edit_counts.py](./analysis/boxb_in_vivo_sequencing/scripts/plot_recorder_edit_counts.py) | [figure_5c.csv](./source_data/figure_5c.csv) |
| 5D | [recorder_sites_and_placements_tada_lambdaN_original](./analysis/boxb_in_vivo_sequencing/figures/recorder_sites_and_placements_tada_lambdaN_original.pdf) | [plot_recorder_site_placement.py](./analysis/boxb_in_vivo_sequencing/scripts/plot_recorder_site_placement.py) | [figure_5d.csv](./source_data/figure_5d.csv) |
| 5E | [boxb_all_variant_A8_vs_A2_original](./analysis/boxb_in_vivo_sequencing/figures/boxb_all_variant_A8_vs_A2_original.pdf) | [plot_boxb_site_correlation.py](./analysis/boxb_in_vivo_sequencing/scripts/plot_boxb_site_correlation.py) | [figure_5e.csv](./source_data/figure_5e.csv) |
| 5F | [boxb_mutation_distance_mean_A8_A2_original](./analysis/boxb_in_vivo_sequencing/figures/boxb_mutation_distance_mean_A8_A2_original.pdf) | [plot_boxb_mutation_distance.py](./analysis/boxb_in_vivo_sequencing/scripts/plot_boxb_mutation_distance.py) | [figure_5f.csv](./source_data/figure_5f.csv) |
| 5G | [boxb_single_mutant_heatmap_mean_A8_A2_original](./analysis/boxb_in_vivo_sequencing/figures/boxb_single_mutant_heatmap_mean_A8_A2_original.pdf) | [plot_boxb_effect_heatmaps.py](./analysis/boxb_in_vivo_sequencing/scripts/plot_boxb_effect_heatmaps.py) | [figure_5g.csv](./source_data/figure_5g.csv) |
| 6B | [variant_read_coverage](./analysis/ln_dms_sequencing/figures/variant_read_coverage.pdf) | [plot_variant_read_coverage.py](./analysis/ln_dms_sequencing/scripts/plot_variant_read_coverage.py) | [figure_6b.csv](./source_data/figure_6b.csv) |
| 6C | [recorder_site_profile](./analysis/ln_dms_sequencing/figures/recorder_site_profile.pdf) | [plot_recorder_site_profile.py](./analysis/ln_dms_sequencing/scripts/plot_recorder_site_profile.py) | [figure_6c.csv](./source_data/figure_6c.csv) |
| 6D | [amino_acid_site_correlation_A2_A15](./analysis/ln_dms_sequencing/figures/amino_acid_site_correlation_A2_A15.pdf) | [plot_amino_acid_site_correlation.py](./analysis/ln_dms_sequencing/scripts/plot_amino_acid_site_correlation.py) | [figure_6d.csv](./source_data/figure_6d.csv) |
| 6E | [variant_class_distribution](./analysis/ln_dms_sequencing/figures/variant_class_distribution.pdf) | [plot_variant_class_distribution.py](./analysis/ln_dms_sequencing/scripts/plot_variant_class_distribution.py) | [figure_6e.csv](./source_data/figure_6e.csv) |
| 6F | [position_editing_profile](./analysis/ln_dms_sequencing/figures/position_editing_profile.pdf) | [plot_position_editing_profile.py](./analysis/ln_dms_sequencing/scripts/plot_position_editing_profile.py) | [figure_6f.csv](./source_data/figure_6f.csv) |
| 6G | [boxb_lambdan_structure_dms_positions](./analysis/boxb_structure_views/figures/boxb_lambdan_structure_dms_positions.pdf) | as above, with the λN positions in [dms_position_summary.csv](./analysis/boxb_structure_views/annotations/dms_position_summary.csv) | — |
| 6H | [amino_acid_editing_logo](./analysis/ln_dms_sequencing/figures/amino_acid_editing_logo.pdf) | [plot_amino_acid_editing_logo.py](./analysis/ln_dms_sequencing/scripts/plot_amino_acid_editing_logo.py) | [figure_6h.csv](./source_data/figure_6h.csv) |
| S1A | [boxb_hairpin_site_editing_enzyme_concentration](./analysis/boxb_in_vitro_sequencing/figures/boxb_hairpin_site_editing_enzyme_concentration.pdf) | [plot_boxb_hairpin_sites.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_hairpin_sites.py) | [figure_s1a.csv](./source_data/figure_s1a.csv) |
| S1B | [recorder_edit_counts_tada_lambdaN_time](./analysis/boxb_in_vitro_sequencing/figures/recorder_edit_counts_tada_lambdaN_time.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_s1b.csv](./source_data/figure_s1b.csv) |
| S2A | [recorder_site_distance_tada_lambdaN_125nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_distance_tada_lambdaN_125nM.pdf) <br> [recorder_site_distance_tada_125nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_distance_tada_125nM.pdf) | [plot_recorder_site_distances.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_site_distances.py) | [figure_s2a.csv](./source_data/figure_s2a.csv) |
| S2B | [recorder_site_distance_tada_lambdaN_500nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_distance_tada_lambdaN_500nM.pdf) <br> [recorder_site_distance_tada_500nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_distance_tada_500nM.pdf) | [plot_recorder_site_distances.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_site_distances.py) | [figure_s2b.csv](./source_data/figure_s2b.csv) |
| S3A | [recorder_context_editing_sites](./analysis/boxb_in_vitro_sequencing/figures/recorder_context_editing_sites.pdf) | [plot_recorder_contexts.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_recorder_contexts.py) | [figure_s3a.csv](./source_data/figure_s3a.csv) |
| S3B | [rabe_rna_sequence_context](./analysis/published_rabe_contexts/figures/rabe_rna_sequence_context.pdf) | [plot_context.r](./analysis/published_rabe_contexts/scripts/plot_context.r) | [figure_s3b.csv](./source_data/figure_s3b.csv) |
| S4A | [boxb_all_variant_A8_vs_A2_125nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_all_variant_A8_vs_A2_125nM.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_s4a.csv](./source_data/figure_s4a.csv) |
| S4B | [boxb_mutation_distance_editing_mean_A8_A2_125nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_mutation_distance_editing_mean_A8_A2_125nM.pdf) | [plot_boxb_mutation_distance.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_mutation_distance.py) | [figure_s4b.csv](./source_data/figure_s4b.csv) |
| S4C | [boxb_single_mutant_editing_mean_A8_A2_125nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_single_mutant_editing_mean_A8_A2_125nM.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_s4c.csv](./source_data/figure_s4c.csv) |
| S4D | [boxb_all_variant_A8_vs_A2_500nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_all_variant_A8_vs_A2_500nM.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_s4d.csv](./source_data/figure_s4d.csv) |
| S4E | [boxb_mutation_distance_editing_mean_A8_A2_500nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_mutation_distance_editing_mean_A8_A2_500nM.pdf) | [plot_boxb_mutation_distance.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_mutation_distance.py) | [figure_s4e.csv](./source_data/figure_s4e.csv) |
| S4F | [boxb_single_mutant_editing_mean_A8_A2_500nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_single_mutant_editing_mean_A8_A2_500nM.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_s4f.csv](./source_data/figure_s4f.csv) |
| S5A | [boxb_multimutant_heatmaps_mean_A8_A2_lambdaN_125nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_multimutant_heatmaps_mean_A8_A2_lambdaN_125nM.pdf) | [plot_boxb_multimutant_heatmaps.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_multimutant_heatmaps.py) | [figure_s5a.csv](./source_data/figure_s5a.csv) |
| S5B | [boxb_multimutant_heatmaps_mean_A8_A2_lambdaN_500nM](./analysis/boxb_in_vitro_sequencing/figures/boxb_multimutant_heatmaps_mean_A8_A2_lambdaN_500nM.pdf) | [plot_boxb_multimutant_heatmaps.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_multimutant_heatmaps.py) | [figure_s5b.csv](./source_data/figure_s5b.csv) |
| S6A | [boxb_single_mutant_editing_mean_A8_A2_time](./analysis/boxb_in_vitro_sequencing/figures/boxb_single_mutant_editing_mean_A8_A2_time.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_s6a.csv](./source_data/figure_s6a.csv) |
| S6B | [boxb_mutation_distance_editing_mean_A8_A2_time](./analysis/boxb_in_vitro_sequencing/figures/boxb_mutation_distance_editing_mean_A8_A2_time.pdf) | [plot_boxb_mutation_distance.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_mutation_distance.py) | [figure_s6b.csv](./source_data/figure_s6b.csv) |
| S6C | [boxb_all_variant_A8_vs_A2_time](./analysis/boxb_in_vitro_sequencing/figures/boxb_all_variant_A8_vs_A2_time.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_s6c.csv](./source_data/figure_s6c.csv) |
| S6D | [boxb_single_mutant_A8_vs_A2_time](./analysis/boxb_in_vitro_sequencing/figures/boxb_single_mutant_A8_vs_A2_time.pdf) | [plot_boxb_substitutions.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_substitutions.py) | [figure_s6d.csv](./source_data/figure_s6d.csv) |
| S7B | [recorder_site_editing_pag_components_A2_250nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_editing_pag_components_A2_250nM.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_s7b.csv](./source_data/figure_s7b.csv) |
| S7C | [recorder_site_editing_editors_A8_250nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_editing_editors_A8_250nM.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_s7c.csv](./source_data/figure_s7c.csv) |
| S7D | [recorder_site_editing_editors_A2_250nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_editing_editors_A2_250nM.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_s7d.csv](./source_data/figure_s7d.csv) |
| S7E | [recorder_site_editing_pag_conditions_A2](./analysis/boxb_in_vitro_sequencing/figures/recorder_site_editing_pag_conditions_A2.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_s7e.csv](./source_data/figure_s7e.csv) |
| S7F | [recorder_edit_counts_pag_conditions](./analysis/boxb_in_vitro_sequencing/figures/recorder_edit_counts_pag_conditions.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_s7f.csv](./source_data/figure_s7f.csv) |
| S7G | [recorder_edit_counts_pag_components_250nM](./analysis/boxb_in_vitro_sequencing/figures/recorder_edit_counts_pag_components_250nM.pdf) | [plot_enzyme_comparisons.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_enzyme_comparisons.py) | [figure_s7g.csv](./source_data/figure_s7g.csv) |
| S7H | [boxb_direct_vs_split_recruitment_mean_A8_A2](./analysis/boxb_in_vitro_sequencing/figures/boxb_direct_vs_split_recruitment_mean_A8_A2.pdf) | [plot_boxb_recruitment_comparison.py](./analysis/boxb_in_vitro_sequencing/scripts/plot_boxb_recruitment_comparison.py) | [figure_s7h.csv](./source_data/figure_s7h.csv) |
| S8A | [recorder_edit_counts_tada_lambdaN_vs_no_editor_original](./analysis/boxb_in_vivo_sequencing/figures/recorder_edit_counts_tada_lambdaN_vs_no_editor_original.pdf) | [plot_recorder_edit_counts.py](./analysis/boxb_in_vivo_sequencing/scripts/plot_recorder_edit_counts.py) | [figure_s8a.csv](./source_data/figure_s8a.csv) |
| S8B | [boxb_in_vitro_vs_in_vivo_mean_A8_A2](./analysis/boxb_in_vivo_sequencing/figures/boxb_in_vitro_vs_in_vivo_mean_A8_A2.pdf) | [plot_in_vitro_comparison.py](./analysis/boxb_in_vivo_sequencing/scripts/plot_in_vitro_comparison.py) | [figure_s8b.csv](./source_data/figure_s8b.csv) |

## Reproduce

The analyses are Snakemake workflows that run in Singularity containers, launched by
[run_everything.sh](./run_everything.sh) from the repository root:

```bash
./run_everything.sh
```

```bash
./run_everything.sh --plots
```

The first command counts reads and then redraws; the second redraws figures from count tables that
are already present. `--only <name>` restricts either to one analysis, and a trailing rule name runs
one figure. Set `SNAKEMAKE`, `SINGULARITY_PREFIX`, `SINGULARITY_ARGS` and `CORES` to suit your
machine.

Raw reads come from the BioProject, fetched by
[download_from_sra.sh](./data/download_from_sra.sh):

```bash
./data/download_from_sra.sh
```

This fetches each run by accession and names the files as the workflows expect. Counting from FASTQ
needs a compute allocation and several hours, and writes about 130 GB into the ignored data
directory of each analysis folder; redrawing takes minutes.

Containers are
[`ghcr.io/rasilab/r_python:2.4.2`](https://github.com/orgs/rasilab/packages/container/package/r_python)
throughout, except the in vitro read counter, which uses
[`ghcr.io/rasilab/python:1.0.0`](https://github.com/orgs/rasilab/packages/container/package/python),
and the published-data reanalysis, which uses
[`ghcr.io/rasilab/hisat2_3n:1.0.1`](https://github.com/orgs/rasilab/packages/container/package/hisat2_3n).

## Layout

| Folder | Contents |
| --- | --- |
| [nanopore_reporter_editing](./analysis/nanopore_reporter_editing/README.md) | Low-throughput reporter editing, Figure 1C |
| [boxb_in_vitro_sequencing](./analysis/boxb_in_vitro_sequencing/README.md) | Purified-enzyme boxB library, Figures 2-4 and Supplementary 1-7 |
| [boxb_in_vivo_sequencing](./analysis/boxb_in_vivo_sequencing/README.md) | Cellular boxB library, Figure 5 and Supplementary 8 |
| [ln_dms_sequencing](./analysis/ln_dms_sequencing/README.md) | λN deep mutational scan, Figure 6 |
| [published_rabe_contexts](./analysis/published_rabe_contexts/README.md) | Reanalysis of published rABE RNA sequencing, Supplementary 3B |
| [boxb_structure_views](./analysis/boxb_structure_views/README.md) | λN-boxB structure panels, Figures 3E and 6G |

Within each analysis folder, annotations holds the versioned CSVs that define samples, read layouts,
coordinates and figure selections; scripts holds the workflows and the plotting code; figures holds
the output; tables holds summary tables; and data holds reads and count tables, which are not
versioned.

| Summary table | Contents |
| --- | --- |
| [figure_1c_statistics.csv](./analysis/nanopore_reporter_editing/tables/figure_1c_statistics.csv) | Fold changes and Fisher exact tests behind Figure 1C |
| [wildtype_statistics.csv](./analysis/ln_dms_sequencing/tables/wildtype_statistics.csv) | λN wild-type dispersion, depth contrast, test calibration and robustness series |
| [rabe_rna_sequence_context.csv](./analysis/published_rabe_contexts/tables/rabe_rna_sequence_context.csv) | Per-context editing in the published rABE samples |
