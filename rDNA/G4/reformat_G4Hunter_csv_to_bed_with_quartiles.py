#!/usr/bin/env python3

import csv
import sys
import numpy as np

if len(sys.argv) != 5:
    sys.exit(
        f"Usage: {sys.argv[0]} input.tsv output.bed min_abs_score sequence_name"
    )

infile = sys.argv[1]
outfile = sys.argv[2]
threshold = float(sys.argv[3])
seqname = sys.argv[4]

# Read and filter
records = []

with open(infile, newline="") as fin:
    reader = csv.DictReader(fin, delimiter="\t")

    for row in reader:
        score = float(row["SCORE"].strip('"'))

        if abs(score) < threshold:
            continue

        start = int(row["POSITION"].strip('"')) - 1
        length = int(row["LENGTH"].strip('"'))
        end = start + length

        records.append({
            "start": start,
            "end": end,
            "score": score
        })

if not records:
    sys.exit(
        f"ERROR: No records remain after filtering abs(SCORE) >= {threshold}"
    )

# Calculate quartiles on filtered absolute scores
abs_scores = np.array([abs(r["score"]) for r in records])

q25, q50, q75 = np.percentile(abs_scores, [25, 50, 75])

def assign_quartile(value):
    if value <= q25:
        return "Q1"
    elif value <= q50:
        return "Q2"
    elif value <= q75:
        return "Q3"
    else:
        return "Q4"

# Write BED
with open(outfile, "w") as fout:
    for r in records:
        strand = "-" if r["score"] < 0 else "+"
        q = assign_quartile(abs(r["score"]))

        fout.write(
            f"{seqname}\t{r['start']}\t{r['end']}\t{q}\t{r['score']}\t{strand}\n"
        )

print(f"Number of retained regions: {len(records)}")
print(f"Absolute score threshold: {threshold}")
print("Quartile boundaries:")
print(f"  Q1/Q2: {q25}")
print(f"  Q2/Q3: {q50}")
print(f"  Q3/Q4: {q75}")