#!/usr/bin/env bash
# Reproduce every figure, source-data CSV and summary table in this repository.
#
# Usage, from the repository root:
#   ./run_everything.sh                    # run what is out of date
#   ./run_everything.sh --plots            # skip read counting and redraw from existing count tables
#   ./run_everything.sh --forceall         # recount and redraw everything from the raw reads
#   ./run_everything.sh --only ln_dms      # just one analysis, by any part of its folder name
#   ./run_everything.sh --plots --only boxb_in_vivo plot_recorder_edit_counts   # one rule
#
# Arguments after the options are passed through to snakemake, so a rule or an output file runs
# just that target.
#
# Counting reads from FASTQ needs a compute allocation and several hours; redrawing takes minutes.
# Count tables for --plots come from ./data/download_count_tables.sh, raw reads from
# ./data/download_from_sra.sh.
#
# Set SNAKEMAKE to the snakemake executable, SINGULARITY_PREFIX to the directory that caches the
# container images, and CORES to the number of cores to use.
set -euo pipefail

SNAKEMAKE="${SNAKEMAKE:-snakemake}"
SINGULARITY_PREFIX="${SINGULARITY_PREFIX:-.singularity}"
SINGULARITY_ARGS="${SINGULARITY_ARGS:-}"
CORES="${CORES:-4}"

COUNTING=(
    analysis/boxb_in_vitro_sequencing/scripts/run_analysis.smk
    analysis/boxb_in_vivo_sequencing/scripts/run_analysis.smk
    analysis/ln_dms_sequencing/scripts/run_counts.smk
    analysis/published_rabe_contexts/scripts/run_analysis.smk
)
PLOTTING=(
    analysis/nanopore_reporter_editing/scripts/plot_figures.smk
    analysis/boxb_in_vitro_sequencing/scripts/plot_figures.smk
    analysis/boxb_in_vivo_sequencing/scripts/plot_figures.smk
    analysis/ln_dms_sequencing/scripts/plot_figures.smk
)
# One count table each plotting workflow reads, as a marker that its tables are already here. The
# Figure 1C workflow has no marker, because its counts are versioned in the repository.
declare -A COUNT_MARKER=(
    [analysis/boxb_in_vitro_sequencing/scripts/plot_figures.smk]=analysis/boxb_in_vitro_sequencing/data/summary_stats_combined/i79_p2.csv.gz
    [analysis/boxb_in_vivo_sequencing/scripts/plot_figures.smk]=analysis/boxb_in_vivo_sequencing/data/recorder_counts/quality/i79_p21.csv.gz
    [analysis/ln_dms_sequencing/scripts/plot_figures.smk]=analysis/ln_dms_sequencing/data/recorder_counts/i84p1_reads.csv
)

if [[ ! -f run_everything.sh ]]; then
    printf 'Run this launcher from the repository root.\n' >&2
    exit 1
fi

workflows=("${COUNTING[@]}" "${PLOTTING[@]}")
plots_only=0
if [[ "${1:-}" == "--plots" ]]; then
    workflows=("${PLOTTING[@]}")
    plots_only=1
    shift
fi
if [[ "${1:-}" == "--only" ]]; then
    [[ -n "${2:-}" ]] || { printf -- '--only needs part of an analysis folder name.\n' >&2; exit 1; }
    selected=()
    for workflow in "${workflows[@]}"; do
        [[ "$workflow" == *"$2"* ]] && selected+=("$workflow")
    done
    [[ ${#selected[@]} -gt 0 ]] || { printf 'No workflow matches %s.\n' "$2" >&2; exit 1; }
    workflows=("${selected[@]}")
    shift 2
fi
extra=("$@")

if [[ $plots_only -eq 1 ]]; then
    for workflow in "${workflows[@]}"; do
        marker="${COUNT_MARKER[$workflow]:-}"
        if [[ -n "$marker" && ! -f "$marker" ]]; then
            printf '%s reads count tables that are not here, starting with %s.\n' "$workflow" "$marker" >&2
            printf 'Fetch them with ./data/download_count_tables.sh, or count them from the raw reads\n' >&2
            printf 'with ./data/download_from_sra.sh and ./run_everything.sh without --plots.\n' >&2
            exit 1
        fi
    done
fi

for workflow in "${workflows[@]}"; do
    printf '\n=== %s\n' "$workflow"
    "$SNAKEMAKE" --snakefile "$workflow" --cores "$CORES" --use-singularity \
        --singularity-prefix "$SINGULARITY_PREFIX" \
        ${SINGULARITY_ARGS:+--singularity-args "$SINGULARITY_ARGS"} \
        --rerun-incomplete "${extra[@]+"${extra[@]}"}"
done
printf '\nAll workflows complete. Figures are in each analysis folder; panel data is in source_data/.\n'
