#!/usr/bin/env python3
import csv
import os
import sys

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <metadata.csv>")
    sys.exit(1)

csv_path = sys.argv[1]

missing_fastq = []
invalid_controls = []
ids = set()
rows = []

# Read CSV and collect IDs
with open(csv_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        ids.add(row['id'])
        rows.append(row)

# Check fastq paths and control validity
for row in rows:
    fq1 = row.get('fastq_path_1', '').strip()
    fq2 = row.get('fastq_path_2', '').strip()
    if not os.path.exists(fq1):
        missing_fastq.append((row['id'], 'fastq_path_1', fq1))
    if not os.path.exists(fq2):
        missing_fastq.append((row['id'], 'fastq_path_2', fq2))
    control = row.get('control', '').strip()
    # Only check control validity if control is not empty
    if control != '' and control not in ids:
        invalid_controls.append((row['id'], control))

# Report
if missing_fastq:
    print("Missing FASTQ files:")
    for sid, col, path in missing_fastq:
        print(f"  Sample {sid}: {col} -> {path} (NOT FOUND)")
else:
    print("All FASTQ files exist.")

if invalid_controls:
    print("Invalid control IDs:")
    for sid, ctrl in invalid_controls:
        print(f"  Sample {sid}: control -> {ctrl} (NOT FOUND in id column)")
else:
    print("All control IDs are valid.")

if not missing_fastq and not invalid_controls:
    print("\nMetadata CSV validation PASSED.")
else:
    print("\nMetadata CSV validation FAILED.")
    sys.exit(2)
