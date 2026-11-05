#!/usr/bin/env bash
# Fetch the count tables that redraw every figure without counting reads from FASTQ.
#
# Usage, from the repository root:
#   ./data/download_count_tables.sh
#
# The archive holds the count tables the four plotting workflows read, and unpacks them into the
# ignored data directory of each analysis folder. Figure 1C needs none of them: its counts are
# versioned in the repository.
#
# Needs curl and tar.
set -euo pipefail

RELEASE="count-tables-v1"
ARCHIVE="count_tables.tar.gz"
URL="https://github.com/rasilab/tada_boxb/releases/download/$RELEASE/$ARCHIVE"
SHA256="f60748305d16a2b07f24ceb53ef1875b176269310e0f0e1b4db65ea7a9ce7e8c"
DESTINATION="data/$ARCHIVE"

[[ -f run_everything.sh ]] || { printf 'Run this from the repository root.\n' >&2; exit 1; }
command -v curl > /dev/null || { printf 'curl is not on PATH.\n' >&2; exit 1; }

if [[ ! -f "$DESTINATION" ]]; then
    printf 'downloading  %s\n' "$URL"
    curl --location --fail --output "$DESTINATION" "$URL"
fi
printf '%s  %s\n' "$SHA256" "$DESTINATION" | sha256sum --check --quiet

tar -xzf "$DESTINATION"
printf '\nUnpacked %s count tables. Now run ./run_everything.sh --plots\n' "$(tar -tzf "$DESTINATION" | wc -l)"
