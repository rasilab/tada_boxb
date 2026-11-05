"""Refresh public run metadata only; never download sequencing reads."""

import csv
import datetime
import io
import os
import re
import urllib.parse
import urllib.request

ANALYSIS = "analysis/published_rabe_contexts"
SAMPLES = f"{ANALYSIS}/annotations/sample_annotations.csv"
OUTPUT = f"{ANALYSIS}/annotations/run_annotations.csv"
ENDPOINT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = ["run_accession", "experiment_accession", "sample_accession", "study_accession",
          "library_layout", "read_count", "base_count", "fastq_ftp", "fastq_md5"]
TIMEOUT = 60

with open(SAMPLES, newline="") as handle:
    samples = list(csv.DictReader(handle))
rows = []
for sample in samples:
    query = urllib.parse.urlencode({"accession": sample["experiment_accession"],
                                    "result": "read_run", "fields": ",".join(FIELDS),
                                    "format": "tsv"})
    url = f"{ENDPOINT}?{query}"
    with urllib.request.urlopen(url, timeout=TIMEOUT) as response:
        records = list(csv.DictReader(io.StringIO(response.read().decode()), delimiter="\t"))
    if not records:
        raise ValueError(f"No public runs for {sample['experiment_accession']}")
    for record in records:
        if record["experiment_accession"] != sample["experiment_accession"]:
            raise ValueError("Public metadata does not match the requested experiment")
        links = record["fastq_ftp"].split(";")
        checksums = record["fastq_md5"].split(";")
        if record["library_layout"] != "PAIRED" or len(links) != 2 or len(checksums) != 2:
            raise ValueError("Unexpected library layout or FASTQ set; inspect metadata before proceeding")
        mates = dict(zip(links, checksums))
        row = {"sample_name": sample["sample_name"], **record}
        for mate in (1, 2):
            matches = [link for link in links if link.endswith(f"_{mate}.fastq.gz")]
            if len(matches) != 1 or not re.fullmatch(r"[0-9a-f]{32}", mates[matches[0]]):
                raise ValueError("Missing or ambiguous FASTQ mate/checksum")
            row[f"read{mate}_url"] = f"https://{matches[0]}"
            row[f"read{mate}_md5"] = mates[matches[0]]
        row["metadata_url"] = url
        row["retrieved_date"] = datetime.date.today().isoformat()
        rows.append(row)
if len({row["run_accession"] for row in rows}) != len(rows):
    raise ValueError("A run was assigned more than once")
with open(f"{OUTPUT}.partial", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
os.replace(f"{OUTPUT}.partial", OUTPUT)
print(f"Recorded metadata for {len(rows)} runs from {len(samples)} samples; no reads downloaded")
