#!/usr/bin/env bash
# Reproduce every figure, source-data CSV and summary table in this repository.
#
# Usage, from the repository root:
#   ./run_everything.sh                    # run what is out of date
#   ./run_everything.sh --plots            # skip read counting and redraw from existing count tables
#   ./run_everything.sh --singularity      # run every step inside the container it declares
#   ./run_everything.sh --only ln_dms      # just one analysis, by any part of its folder name
#   ./run_everything.sh --forceall         # recount and redraw everything from the raw reads
#   ./run_everything.sh --plots --only boxb_in_vivo plot_recorder_edit_counts   # one rule
#
# The options above may be given in any order. Anything after them is passed through to snakemake,
# so a rule name or an output file runs just that target.
#
# By default every step runs with whatever is on PATH, which needs only a Python environment holding
# pandas, numpy, matplotlib, scipy and biopython. --singularity instead runs each step inside the
# container its workflow declares, which is how the published figures were produced, and is the
# easiest way to run the published-data reanalysis, whose aligner is otherwise yours to install.
#
# Counting reads from FASTQ needs a compute allocation and several hours; redrawing takes minutes.
# Count tables for --plots come from ./data/download_count_tables.sh, raw reads from
# ./data/download_from_sra.sh.
#
# Set SNAKEMAKE to the snakemake executable and CORES to the number of cores to use. With
# --singularity, set SINGULARITY_PREFIX to the directory that caches the container images and
# SINGULARITY_ARGS to the bind mounts this machine needs. The README tabulates all four.
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
# Programs a workflow calls out to when it runs without containers. Python itself is needed
# everywhere and comes with snakemake; only the published-data reanalysis runs separate binaries.
declare -A REQUIRED_PROGRAMS=(
    [analysis/published_rabe_contexts/scripts/run_analysis.smk]="hisat-3n hisat-3n-build hisat-3n-table samtools"
)

if [[ ! -f run_everything.sh ]]; then
    printf 'Run this launcher from the repository root.\n' >&2
    exit 1
fi

workflows=("${COUNTING[@]}" "${PLOTTING[@]}")
plots_only=0
use_singularity=0
only=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --plots) plots_only=1; shift ;;
        --singularity) use_singularity=1; shift ;;
        --only)
            [[ -n "${2:-}" ]] || { printf -- '--only needs part of an analysis folder name.\n' >&2; exit 1; }
            only="$2"
            shift 2 ;;
        *) break ;;
    esac
done
if [[ $plots_only -eq 1 ]]; then
    workflows=("${PLOTTING[@]}")
fi
if [[ -n "$only" ]]; then
    selected=()
    for workflow in "${workflows[@]}"; do
        [[ "$workflow" == *"$only"* ]] && selected+=("$workflow")
    done
    [[ ${#selected[@]} -gt 0 ]] || { printf 'No workflow matches %s.\n' "$only" >&2; exit 1; }
    workflows=("${selected[@]}")
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

singularity_flags=()
if [[ $use_singularity -eq 1 ]]; then
    singularity_flags=(--use-singularity --singularity-prefix "$SINGULARITY_PREFIX")
    if [[ -n "$SINGULARITY_ARGS" ]]; then
        singularity_flags+=(--singularity-args "$SINGULARITY_ARGS")
    fi
else
    for workflow in "${workflows[@]}"; do
        for program in ${REQUIRED_PROGRAMS[$workflow]:-}; do
            if ! command -v "$program" > /dev/null; then
                printf '%s runs %s, which is not on PATH.\n' "$workflow" "$program" >&2
                printf 'Install it, or rerun with --singularity to take it from the container instead.\n' >&2
                exit 1
            fi
        done
    done
fi

for workflow in "${workflows[@]}"; do
    printf '\n=== %s\n' "$workflow"
    "$SNAKEMAKE" --snakefile "$workflow" --cores "$CORES" \
        ${singularity_flags[@]+"${singularity_flags[@]}"} \
        --rerun-incomplete "${extra[@]+"${extra[@]}"}"
done
printf '\nAll workflows complete. Figures are in each analysis folder; panel data is in source_data/.\n'
