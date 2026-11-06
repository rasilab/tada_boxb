#!/usr/bin/env bash
# Fetch the raw reads from the BioProject and lay them out where the workflows expect them.
#
# Usage, from the repository root:
#   ./data/download_from_sra.sh                # every run, with fasterq-dump from PATH
#   ./data/download_from_sra.sh i79p1 i79p3    # only these, by the sample_id column
#   ./data/download_from_sra.sh --singularity  # run fasterq-dump inside its container instead
#
# data/sra_accessions.csv lists the runs in BioProject PRJNA1255650 against the sample_id each
# analysis uses. Its read1_fastq column is where read 1 lands, named as the workflows expect;
# read 2 of a paired run swaps _R1_ for _R2_.
#
# By default this calls fasterq-dump from PATH, which conda installs with
# `conda install -c bioconda sra-tools`. --singularity instead runs it inside the container named
# below, so that nothing but Singularity has to be installed: set SINGULARITY to apptainer if that
# is what you have, SINGULARITY_PREFIX to the directory that caches container images, and
# SINGULARITY_ARGS to the bind mounts your machine needs. THREADS sets how many threads each
# download may use. These are the same knobs run_everything.sh takes.
set -euo pipefail

ACCESSIONS="data/sra_accessions.csv"
CONTAINER="docker://ghcr.io/rasilab/sratools:3.0.8"
IMAGE_NAME="sratools_3.0.8.sif"
SINGULARITY="${SINGULARITY:-singularity}"
SINGULARITY_PREFIX="${SINGULARITY_PREFIX:-.singularity}"
SINGULARITY_ARGS="${SINGULARITY_ARGS:-}"
THREADS="${THREADS:-4}"

use_singularity=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --singularity) use_singularity=1; shift ;;
        *) break ;;
    esac
done

[[ -f "$ACCESSIONS" ]] || { printf 'Run this from the repository root.\n' >&2; exit 1; }

if [[ $use_singularity -eq 1 ]]; then
    command -v "$SINGULARITY" > /dev/null || {
        printf '%s is not on PATH. Install Singularity or Apptainer, or drop --singularity and\n' "$SINGULARITY" >&2
        printf 'install sra-tools instead.\n' >&2
        exit 1; }
    IMAGE="$SINGULARITY_PREFIX/$IMAGE_NAME"
    mkdir -p "$SINGULARITY_PREFIX"
    if [[ ! -e "$IMAGE" ]]; then
        printf 'building         %s\n' "$IMAGE"
        "$SINGULARITY" build "$IMAGE" "$CONTAINER"
    fi
    fasterq=("$SINGULARITY" exec ${SINGULARITY_ARGS:+$SINGULARITY_ARGS} "$IMAGE" fasterq-dump)
else
    command -v fasterq-dump > /dev/null || {
        printf 'fasterq-dump is not on PATH. Install sra-tools with\n' >&2
        printf '    conda install -c bioconda sra-tools\n' >&2
        printf 'or rerun with --singularity to take it from a container instead.\n' >&2
        exit 1; }
    fasterq=(fasterq-dump)
fi

wanted=("$@")
while IFS=, read -r analysis sample run _experiment _biosample layout read1; do
    [[ "$analysis" == \#* || "$analysis" == "analysis" || -z "$analysis" ]] && continue
    if [[ ${#wanted[@]} -gt 0 ]] && ! printf '%s\n' "${wanted[@]}" | grep -qx "$sample"; then
        continue
    fi
    if [[ -f "$read1" ]]; then
        printf 'already present  %s\n' "$read1"
        continue
    fi
    directory="$(dirname "$read1")"
    unpacked="${read1%.gz}"
    mkdir -p "$directory"
    printf 'downloading      %s  %s\n' "$sample" "$run"
    "${fasterq[@]}" "$run" --split-files --include-technical --threads "$THREADS" \
        --outdir "$directory"
    # fasterq-dump numbers paired streams in submission order, R1 then R2, and names a single
    # stream after the run alone.
    if [[ "$layout" == "paired" ]]; then
        mv "$directory/${run}_1.fastq" "$unpacked"
        mv "$directory/${run}_2.fastq" "${unpacked/_R1_/_R2_}"
    else
        mv "$directory/${run}.fastq" "$unpacked"
    fi
    if [[ "$read1" == *.gz ]]; then
        gzip -f "$unpacked"
    fi
done < "$ACCESSIONS"

printf '\nDone. Now run ./run_everything.sh\n'
