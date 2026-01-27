#!/usr/bin/env python3
# usage:
#   python classify_midpoints_dbscan.py [input.bed] [eps]
# examples:
#   python classify_midpoints_dbscan.py chr1.bed
#   python classify_midpoints_dbscan.py chr1.bed 300000
#   cat chr1.bed | python classify_midpoints_dbscan.py 300000
#
# This script classifies BED intervals using DBSCAN on midpoints.
# Steps:
#   1) Compute interval lengths.
#   2) Identify "long intervals": top half by length AND length >= 10000.
#   3) Compute mean & std dev of midpoints of those long intervals.
#   4) Compute z-scores for ALL midpoints relative to that mean/std.
#   5) Points >3σ are outliers (smaller DBSCAN weights).
#   6) Use weighted DBSCAN on all midpoints.
#
# Only rows with dbscan_label != -1 are written to stdout, the outliers go to stderr
#
# cechova.biomonika@gmail.com

import sys
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

# --- Arguments ---
if len(sys.argv) == 1:
    bed_path = None
    eps = 300000.0
elif len(sys.argv) == 2:
    arg = sys.argv[1]
    try:
        eps = float(arg)
        bed_path = None
    except ValueError:
        bed_path = arg
        eps = 300000.0
else:
    bed_path = sys.argv[1]
    eps = float(sys.argv[2])

min_samples = 1  # default

# --- Read BED ---
if bed_path:
    df = pd.read_csv(bed_path, sep="\t", header=None, comment="#")
else:
    if sys.stdin.isatty():
        sys.stderr.write("Error: No input provided (expected a BED file or piped data)\n")
        sys.exit(1)
    df = pd.read_csv(sys.stdin, sep="\t", header=None, comment="#")

if df.shape[1] < 3:
    sys.stderr.write("BED must have at least 3 columns: chrom, start, end\n")
    sys.exit(1)

# Ensure numeric
df[1] = df[1].astype(int)
df[2] = df[2].astype(int)

# --- Compute lengths & midpoints ---
lengths = (df[2] - df[1]).to_numpy()
midpoints = ((df[1] + df[2]) // 2).to_numpy().astype(float)

if len(df) == 0:
    sys.stderr.write("No intervals found.\n")
    sys.exit(0)

# --- Identify long intervals ---
n = len(df)
half_n = max(1, n // 2)

# Get indices of top half by length
top_idx = np.argpartition(lengths, -half_n)[-half_n:]
long_mask = np.zeros(n, dtype=bool)
long_mask[top_idx] = True

# Also require minimum length >= 10000
long_mask &= (lengths >= 10000)

# Extract long intervals
long_midpoints = midpoints[long_mask]

# --- Fallback if no intervals meet both criteria ---
if len(long_midpoints) == 0:
    sys.stderr.write("Warning: No intervals meet 'long' criteria (top half & >=10000 bp). Using all intervals instead.\n")
    long_midpoints = midpoints

# --- Compute mean & std from long intervals ---
mean_long = np.mean(long_midpoints)
std_long = np.std(long_midpoints, ddof=0)
if std_long == 0:
    std_long = 100000.0  # avoid divide-by-zero
sys.stderr.write(f"std_long = {std_long}\n")

# --- Compute z-scores for ALL midpoints based on long ones ---
zscores = (midpoints - mean_long) / std_long
is_outlier = np.abs(zscores) > 3
weights = np.where(is_outlier, 0.1, 1.0)

# --- Run DBSCAN ---
X = midpoints.reshape(-1, 1)
db = DBSCAN(eps=eps, min_samples=min_samples)
labels = db.fit_predict(X, sample_weight=weights)

# --- Prepare output ---
out = df.copy()
out[4] = labels  # cluster label

# --- Keep ONLY the cluster with the most elements ---
cluster_labels = out[4][out[4] != -1]
if cluster_labels.empty:
    # No clusters found; send everything to stderr
    sys.stderr.write("Warning: DBSCAN found no clusters (all noise). Outputting all rows to stderr.\n")
    out.to_csv(sys.stderr, sep="\t", header=False, index=False)
    sys.exit(0)

# Find the largest cluster label (by count)
counts = cluster_labels.value_counts()
largest_label = counts.idxmax()
largest_size = counts.loc[largest_label]
sys.stderr.write(f"Largest cluster label = {largest_label} (n={largest_size})\n")

# Split rows
keep = out[out[4] == largest_label]
drop = out[out[4] != largest_label]  # includes other clusters and noise (-1)

# --- Write TSV outputs ---
keep.to_csv(sys.stdout, sep="\t", header=False, index=False)
drop.to_csv(sys.stderr, sep="\t", header=False, index=False)






