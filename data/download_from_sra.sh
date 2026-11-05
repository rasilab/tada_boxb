#!/usr/bin/env bash
# Fetch the raw reads from the BioProject and lay them out where the workflows expect them.
#
# Usage, from the repository root:
#   ./data/download_from_sra.sh              # every run
#   ./data/download_from_sra.sh i79p1 i79p3  # only these, by the sample_id column
#
# data/sra_accessions.csv lists the runs in BioProject PRJNA1255650 against the sample_id each
# analysis uses. Its read1_fastq column is where read 1 lands, named as the workflows expect;
# read 2 of a paired run swaps _R1_ for _R2_.
#
# Needs sra-tools, with fasterq-dump on PATH.
set -euo pipefail

ACCESSIONS="data/sra_accessions.csv"
THREADS="${THREADS:-4}"

[[ -f "$ACCESSIONS" ]] || { printf 'Run this from the repository root.\n' >&2; exit 1; }
command -v fasterq-dump > /dev/null || {
    printf 'fasterq-dump is not on PATH. Install sra-tools first.\n' >&2; exit 1; }

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
    fasterq-dump "$run" --split-files --include-technical --threads "$THREADS" --outdir "$directory"
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
