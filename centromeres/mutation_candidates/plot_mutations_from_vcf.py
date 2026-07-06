import argparse
import sys
import pysam
import matplotlib.pyplot as plt
from collections import Counter
import os
import numpy as np
import matplotlib.patches as mpatches  # <-- needed for custom legend patches

# 12 directional single-nucleotide classes (order matches the background vector)
LABELS = ["A>C","A>G","A>T",
          "C>A","C>G","C>T",
          "G>A","G>C","G>T",
          "T>A","T>C","T>G"]

# Provided background counts (gray bars), in the exact order of LABELS
BACKGROUND_COUNTS = [106, 309, 123, 109, 88, 335, 296, 103, 99, 136, 325, 97]

# Strand-collapsed canonical labels (pyrimidine-centered)
CANONICAL_LABELS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
COMPLEMENT = {"A":"T","C":"G","G":"C","T":"A"}

def to_canonical(mutation):
    ref, alt = mutation.split(">")
    if ref in ("A", "G"):
        ref = COMPLEMENT[ref]
        alt = COMPLEMENT[alt]
    return f"{ref}>{alt}"

def collapse_counts(counts_dict, labels_order):
    collapsed = Counter()
    for k, v in counts_dict.items():
        canon = to_canonical(k)
        if canon in CANONICAL_LABELS:
            collapsed[canon] += v
    return [collapsed.get(k, 0) for k in CANONICAL_LABELS]

def main():
    parser = argparse.ArgumentParser(
        description="Visualize SNV mutation types from a VCF file with background (directional and strand-collapsed)."
    )
    parser.add_argument("--vcf", required=True, help="Path to the input VCF file (.vcf or .vcf.gz)")
    args = parser.parse_args()

    if not os.path.isfile(args.vcf):
        print(f"Error: File '{args.vcf}' does not exist.")
        sys.exit(1)

    try:
        vcf = pysam.VariantFile(args.vcf)
    except Exception as e:
        print(f"Error opening VCF: {e}")
        sys.exit(1)

    if len(BACKGROUND_COUNTS) != len(LABELS):
        print("Error: BACKGROUND_COUNTS length must match number of LABELS (12).")
        sys.exit(1)

    vcf_basename = os.path.basename(args.vcf)
    vcf_prefix = os.path.splitext(vcf_basename)[0]

    # --- Count directional classes from VCF ---
    directional_counts = Counter()
    try:
        iterator = vcf.fetch()
    except Exception:
        iterator = vcf

    for record in iterator:
        ref = record.ref
        alts = record.alts
        if not alts or len(ref) != 1 or ref not in "ACGT":
            continue
        for alt in alts:
            if len(alt) == 1 and alt in "ACGT" and alt != ref:
                m = f"{ref}>{alt}"
                if m in LABELS:
                    directional_counts[m] += 1

    data_values = [directional_counts.get(k, 0) for k in LABELS]
    bg_values = BACKGROUND_COUNTS

    if sum(data_values) == 0:
        print("No SNVs found among the 12 single-nucleotide classes.")
        sys.exit(0)

    # Save observed directional counts (12 labels)
    counts_tsv = f"{vcf_prefix}_observed_directional_counts.tsv"
    with open(counts_tsv, "w") as w:
        w.write("label\tcount\n")
        for lab, cnt in zip(LABELS, data_values):
            w.write(f"{lab}\t{cnt}\n")
    print(f"Observed directional counts saved to: {counts_tsv}")

    # ----------------
    # Plot 1: Counts (Directional)
    # ----------------
    x = np.arange(len(LABELS))
    width = 0.42
    plt.figure(figsize=(12, 6))
    plt.bar(x - width/2, bg_values, width, color="lightgray", edgecolor="black", label="Background")
    plt.bar(x + width/2, data_values, width, color="#990000", edgecolor="black", label="Centromeres")
    plt.xticks(x, LABELS, rotation=45, ha="right", fontsize=14, fontweight="bold")
    plt.xlabel("DNMs: (directional)", fontsize=16)
    plt.ylabel("Count", fontsize=16)
    plt.title("DNMs: Centromere (red) vs Background (gray) [Counts]", fontsize=18)
    plt.legend(fontsize=14)
    plt.tight_layout()
    counts_filename = f"{vcf_prefix}_mutation_counts.png"
    plt.savefig(counts_filename, dpi=300)
    print(f"Counts plot saved to: {counts_filename}")

    # ----------------
    # Plot 2: Percentages (Directional)
    # ----------------
    total_data = sum(data_values)
    total_bg = sum(bg_values)
    data_perc = [(v / total_data * 100) if total_data > 0 else 0 for v in data_values]
    bg_perc = [(v / total_bg * 100) if total_bg > 0 else 0 for v in bg_values]

    plt.figure(figsize=(12, 6))
    plt.bar(x - width/2, bg_perc, width, color="lightgray", edgecolor="black", label="Background")
    plt.bar(x + width/2, data_perc, width, color="#990000", edgecolor="black", label="Centromeres")
    plt.xticks(x, LABELS, rotation=45, ha="right", fontsize=14, fontweight="bold")
    plt.xlabel("DNMs: (directional)", fontsize=16)
    plt.ylabel("Percentage of total (%)", fontsize=16)
    plt.title("DNMs: Centromere (red) vs Background (gray)", fontsize=18)
    plt.legend(fontsize=14)
    plt.tight_layout()
    perc_filename = f"{vcf_prefix}_mutation_percentages.png"
    plt.savefig(perc_filename, dpi=300)
    print(f"Percentages plot saved to: {perc_filename}")

    # ----------------
    # Collapse to 6 canonical classes
    # ----------------
    data_dir_dict = {k: directional_counts.get(k, 0) for k in LABELS}
    data_collapsed = collapse_counts(data_dir_dict, CANONICAL_LABELS)
    bg_dir_dict = {k: v for k, v in zip(LABELS, BACKGROUND_COUNTS)}
    bg_collapsed = collapse_counts(bg_dir_dict, CANONICAL_LABELS)

    # ----------------
    # Plot 3: Counts (Collapsed)
    # ----------------
    x2 = np.arange(len(CANONICAL_LABELS))
    width2 = 0.42
    plt.figure(figsize=(10, 6))
    plt.bar(x2 - width2/2, bg_collapsed, width2, color="lightgray", edgecolor="black", label="Background")
    plt.bar(x2 + width2/2, data_collapsed, width2, color="#990000", edgecolor="black", label="Centromeres")
    plt.xticks(x2, CANONICAL_LABELS, fontsize=14, fontweight="bold")
    plt.xlabel("DNMs: (strand-collapsed)", fontsize=16)
    plt.ylabel("Count", fontsize=16)
    plt.title("DNMs (strand-collapsed): Centromere (red) vs Background (gray)", fontsize=18)
    plt.legend(fontsize=14)
    plt.tight_layout()
    collapsed_counts_filename = f"{vcf_prefix}_mutation_collapsed_counts.png"
    plt.savefig(collapsed_counts_filename, dpi=300)
    print(f"Collapsed counts plot saved to: {collapsed_counts_filename}")

    # ----------------
    # Plot 4: Percentages (Collapsed)
    # ----------------
    total_data_c = sum(data_collapsed)
    total_bg_c = sum(bg_collapsed)
    data_collapsed_perc = [(v / total_data_c * 100) if total_data_c > 0 else 0 for v in data_collapsed]
    bg_collapsed_perc = [(v / total_bg_c * 100) if total_bg_c > 0 else 0 for v in bg_collapsed]

    # define custom 2-item legend (Centromeres & Background)
    centro_patch = mpatches.Patch(color="#990000", label="Centromeres")
    bg_patch     = mpatches.Patch(color="lightgray",  label="Background")

    plt.figure(figsize=(10, 6))
    plt.bar(x2 - width2/2, bg_collapsed_perc, width2, color="lightgray", edgecolor="black")
    plt.bar(x2 + width2/2, data_collapsed_perc, width2, color="#990000", edgecolor="black")
    plt.xticks(x2, CANONICAL_LABELS, fontsize=14, fontweight="bold")
    plt.xlabel("DNMs: (strand-collapsed)", fontsize=16)
    plt.ylabel("Percentage of total (%)", fontsize=16)
    plt.title("DNMs (strand-collapsed): Centromere (red) vs Background (gray)", fontsize=18)
    plt.legend(handles=[centro_patch, bg_patch], fontsize=12, ncol=2)
    plt.tight_layout()
    collapsed_perc_filename = f"{vcf_prefix}_mutation_collapsed_percentages.png"
    plt.savefig(collapsed_perc_filename, dpi=300)
    print(f"Collapsed percentages plot saved to: {collapsed_perc_filename}")

if __name__ == "__main__":
    main()
