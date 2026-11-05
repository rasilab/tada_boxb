"""Combine declared per-barcode count CSVs, preserving the existing count interface."""

import os
import sys
import pandas as pd

OUTPUT_FILE = sys.argv[1]
SAMPLE_ID = sys.argv[2]
INPUT_FILES = sys.argv[3:]
COMPRESSION = {"method": "gzip", "mtime": 0}

frames = []
for filename in sorted(INPUT_FILES):
    data = pd.read_csv(filename)
    data.insert(0, "barcode", os.path.basename(filename).removesuffix(".csv"))
    data.insert(0, "sample_id", SAMPLE_ID)
    frames.append(data)
summary_data = pd.concat(frames, ignore_index=True)
if summary_data.duplicated(["sample_id", "barcode", "insert"]).any():
    raise ValueError("Duplicate count constructs")
summary_data.to_csv(OUTPUT_FILE, index=False, compression=COMPRESSION)
