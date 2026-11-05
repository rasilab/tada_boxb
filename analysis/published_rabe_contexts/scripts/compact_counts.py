"""Stream HISAT-3N's seven-column table into a compressed count CSV."""

import csv
import gzip
import sys

OUTPUT = sys.argv[1]
GZIP_LEVEL = 1
COLUMNS = ["ref", "pos", "strand", "convertedBaseCount", "unconvertedBaseCount"]
EXPECTED = ["ref", "pos", "strand", "convertedBaseQualities", "convertedBaseCount",
            "unconvertedBaseQualities", "unconvertedBaseCount"]

csv.field_size_limit(sys.maxsize)
reader = csv.reader(sys.stdin, delimiter="\t")
if next(reader, None) != EXPECTED:
    raise ValueError("Unexpected HISAT-3N table header")
indices = [EXPECTED.index(column) for column in COLUMNS]
with gzip.open(OUTPUT, "wt", newline="", compresslevel=GZIP_LEVEL) as handle:
    writer = csv.writer(handle)
    writer.writerow(COLUMNS)
    for row in reader:
        if len(row) != len(EXPECTED):
            raise ValueError("Malformed HISAT-3N table row")
        writer.writerow([row[index] for index in indices])
