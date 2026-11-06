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

Everything runs from the repository root through [run_everything.sh](./run_everything.sh), which
launches each analysis as a Snakemake workflow:

```bash
./run_everything.sh
```

```bash
./run_everything.sh --plots
```

The first command counts reads and then redraws; the second redraws figures from count tables that
are already present. `--only <name>` restricts either to one analysis, `--singularity` runs every
step in a container, and a trailing rule name runs one figure. The options may be given in any
order, and anything after them is passed through to Snakemake.

### What has to be installed

Snakemake, and a Python environment holding pandas, numpy, matplotlib, scipy and biopython. By
default every step runs with whatever is on PATH, so nothing further is needed on an ordinary
machine. The published figures were produced under Snakemake 7.0.0 with Python 3.12.8, pandas 2.3.1,
numpy 1.26.4, matplotlib 3.10.3, scipy 1.16.0 and biopython 1.85.

Adding `--singularity` instead runs each step inside the container its workflow declares, so that
the versions above are fixed for you:

```bash
./run_everything.sh --singularity
```

Snakemake then pulls each `docker://` image once into the directory named by `SINGULARITY_PREFIX`
and reuses it from then on, so the first run of an analysis spends a few extra minutes building
images. Apptainer works in place of Singularity.

| Container | Used by |
| --- | --- |
| [`r_python:2.4.2`](https://github.com/orgs/rasilab/packages/container/package/r_python) | every workflow, apart from the two below |
| [`python:1.0.0`](https://github.com/orgs/rasilab/packages/container/package/python) | the in vitro read counter |
| [`hisat2_3n:1.0.1`](https://github.com/orgs/rasilab/packages/container/package/hisat2_3n) | the published-data reanalysis behind Supplementary Figure 3B |
| [`sratools:3.0.8`](https://github.com/orgs/rasilab/packages/container/package/sratools) | [download_from_sra.sh](./data/download_from_sra.sh), which runs `fasterq-dump` inside it |

Every analysis behind the main figures needs nothing but Python. The published-data reanalysis
behind Supplementary Figure 3B is the exception, because it runs HISAT-3N and samtools: give it
`--singularity`, or put those programs on PATH. The launcher checks for them and says which is
missing before starting anything.

Both scripts read the same settings from the environment, and each one has a working default:

| Variable | Read by | Default | Meaning |
| --- | --- | --- | --- |
| `SNAKEMAKE` | [run_everything.sh](./run_everything.sh) | `snakemake` | the Snakemake executable |
| `CORES` | [run_everything.sh](./run_everything.sh) | `4` | cores the workflows may use |
| `THREADS` | [download_from_sra.sh](./data/download_from_sra.sh) | `4` | threads each download may use |
| `SINGULARITY_PREFIX` | both scripts, with `--singularity` | `.singularity` | directory that caches the container images |
| `SINGULARITY_ARGS` | both scripts, with `--singularity` | none | arguments passed through to Singularity, for whatever bind mounts the machine needs |
| `SINGULARITY` | [download_from_sra.sh](./data/download_from_sra.sh), with `--singularity` | `singularity` | set this to `apptainer` instead where that is what is installed |

On a cluster, containers avoid installing anything: bind the filesystem that holds the data, and
point the image cache at shared storage so that the images are built once for everyone.

```bash
export SINGULARITY_PREFIX=/path/to/shared/singularity
export SINGULARITY_ARGS="-B /fh"
export CORES=12
./run_everything.sh --singularity
```

### Count tables, to redraw the figures

Count tables for `--plots` are attached to the
[count-tables-v1](https://github.com/rasilab/tada_boxb/releases/tag/count-tables-v1) release,
fetched by [download_count_tables.sh](./data/download_count_tables.sh):

```bash
./data/download_count_tables.sh
```

That is a 67 MB download, and covers every figure except 1C, whose counts are versioned in the
repository. Redrawing from the tables takes minutes.

### Raw reads, to count from scratch

Raw reads come from the BioProject, fetched by
[download_from_sra.sh](./data/download_from_sra.sh):

```bash
./data/download_from_sra.sh
```

This fetches each run by accession and names the files as the workflows expect. Naming sample
identifiers from the `sample_id` column of [sra_accessions.csv](./data/sra_accessions.csv) fetches
only those runs. It calls `fasterq-dump` from PATH, which conda installs with
`conda install -c bioconda sra-tools`; `--singularity` runs it inside the container instead, so that
nothing has to be installed for it either.

Counting from the reads needs a compute allocation, and far more disk than it downloads:

- The 25 runs are a 6.4 GB download, 20.9 billion bases in all.
- They unpack to 95 GB of FASTQ, because the counting workflows read uncompressed reads. Only the
  λN DMS mRNA run stays gzipped.
- Counting writes per-barcode intermediates on top of that, leaving 167 GB in the ignored data
  directory of each analysis folder, 156 GB of it under the in vitro analysis.
- A complete run of the four analyses behind Figures 1C and 2 to 6 took 1 hour 55 minutes on 12
  cores, divided about evenly between downloading and counting, and never passed 2 GB of memory.

The published-data reanalysis behind Supplementary Figure 3B is separate, and much heavier: building
the HISAT-3N index for GRCh38 alone asks for 300 GB of memory and 16 cores.

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
