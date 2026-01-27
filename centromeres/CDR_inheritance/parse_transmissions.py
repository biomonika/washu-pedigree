import pandas as pd
import os
import glob
import shutil
import subprocess

def process_centromere_transmissions(csv_path):
    """
    Processes a centromere transmissions CSV file and returns filtered inheritance lists.

    Parameters:
        csv_path (str): Path to the centromere transmissions CSV file.

    Returns:
        dict: A dictionary containing two lists:
              - 'two_generational_inheritances'
              - 'three_generational_inheritances'
    """
    # Load the CSV
    centromere_transmissions = pd.read_csv(csv_path, sep=",", header=None)

    # Extract two-generational and three-generational inheritances
    two_generational_inheritances = (
        centromere_transmissions.iloc[1, 1:].dropna().tolist() +
        centromere_transmissions.iloc[2, 1:].dropna().tolist()
    )

    three_generational_inheritances = (
        centromere_transmissions.iloc[9, 1:].dropna().tolist() +
        centromere_transmissions.iloc[11, 1:].dropna().tolist() +
        centromere_transmissions.iloc[13, 1:].dropna().tolist()
    )

    # Handle PAN027-specific chromosomes
    PAN027_three_gen_chroms = [x for x in three_generational_inheritances if "PAN027" in x]
    PAN027_two_gen_chroms = [
        x.replace("maternal", "__TMP__")
         .replace("paternal", "maternal")
         .replace("__TMP__", "paternal")
        for x in PAN027_three_gen_chroms
    ]

    # Filter out PAN027 from the three-generational list
    three_generational_inheritances = [x for x in three_generational_inheritances if "PAN027" not in x]

    # Remove overlaps
    two_generational_inheritances = [
        x for x in two_generational_inheritances if x not in three_generational_inheritances
    ]

    return {
        "two_generational_inheritances": two_generational_inheritances + PAN027_two_gen_chroms,
        "three_generational_inheritances": three_generational_inheritances + PAN027_three_gen_chroms
    }


def loop_by_chromosome(results, src_folder, verbose=True):
    """
    Loop through all chromosomes (chr1..chr22, chrX) and print inheritance info.
    Three-generational lists are ordered as:
        PAN010/PAN011 → PAN027 → PAN028
    """
    two_gen = results["two_generational_inheritances"]
    three_gen = results["three_generational_inheritances"]

    # Chromosomes 1–22 and X
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

    summary = {}

    for chrom in chromosomes:
        pattern = f".{chrom}."

        two = [x for x in two_gen if pattern in x]
        three = [x for x in three_gen if pattern in x]

        # Order three-generational: PAN010/011 → PAN027 → PAN028
        three_sorted = sorted(
            three,
            key=lambda x: (
                0 if "PAN010" in x or "PAN011" in x else
                1 if "PAN027" in x else
                2 if "PAN028" in x else
                3
            )
        )

        summary[chrom] = {"two_gen": two, "three_gen": three_sorted}

        if verbose:
            print(f"\n=== {chrom} ===")
            print(f"Two-generational ({len(two)}): {two}")
            print(f"Three-generational ({len(three_sorted)}): {three_sorted}")

        pan_grandparent_entry = next((x for x in three_sorted if "PAN010" in x or "PAN011" in x), None)
        pan027_entry = next((x for x in three_sorted if "PAN027" in x), None)
        pan028_entry = next((x for x in three_sorted if "PAN028" in x), None)

        # --- Step 1: Copy only relevant files for this chromosome ---
        bed_candidates = glob.glob(os.path.join(src_folder, f"*{chrom}*{pan027_entry}*centrodip.bed"))
        copied_beds = []

        for src in bed_candidates:
            base = os.path.basename(src)
            prefix = base.split(".realigned_to.")[0]  # e.g. "PAN028.chr17.haplotype1"
            if f".{chrom}." in prefix:
                dst = os.path.join(os.getcwd(), base)
                #print("src")
                #print(src)
                #print("dst")
                #print(dst)
                shutil.copy(src, dst)
                copied_beds.append(dst)

        if not copied_beds:
            print(f"  ⚠️ No BED files copied for {chrom}")
            continue

        # --- Step 2: Convert each copied BED file ---
        for bed in copied_beds:
            base = os.path.basename(bed)
            try:
                sample = base.split(".chr")[0]                        # e.g., "PAN028"
                chrom = base.split(".")[1]                            # e.g., "chr17"
                haplotype = base.split(".realigned_to.")[0].split(".")[-1]  # e.g., "haplotype1"
                entity = f"{sample}.{chrom}.{haplotype}"              # e.g., "PAN028.chr17.haplotype1"
            except Exception:
                print(f"  ⚠️ Skipping malformed file: {base}")
                continue

            out_bed = f"{entity}.bed"
            cmd = (
                f'egrep --color "{chrom}" "{base}" | '
                f'grep -v "transition_CDR" | '
                f'awk "{{ if ((\\$3 - \\$2) >= 3000) print }}" | '
                f'python keep_main_dbscan_cluster_filter_outliers.py > "{out_bed}"'
            )
            print(f"  Running: {cmd}")
            subprocess.run(cmd, shell=True, check=False)

        # --- Step 3: Compare three generations (PAN010, PAN027, PAN028) ---
        pan_grandparent = f"{pan_grandparent_entry}.bed"
        pan027 = f"{pan027_entry}.bed"
        pan028 = f"{pan028_entry}.bed"
        if all(os.path.exists(f) for f in [pan_grandparent, pan027, pan028]):
            cmd = (
                'source /opt/miniconda/etc/profile.d/conda.sh && '
                'conda activate /private/groups/migalab/mcechova/conda/genomics-r && '
                f'Rscript compare_gens_overlap_span.R "{pan_grandparent}" "{pan027}" "{pan028}"'
            )

            print(f"  Comparing: {cmd}")

            # Save Rscript output to file named by chromosome (e.g. chr17_results.txt)
            result_file = f"{chrom}_CDR_3G_stats.minL3k.txt"

            with open(result_file, "w") as f:
                subprocess.run(cmd, shell=True, executable="/bin/bash", check=False, stdout=f, stderr=subprocess.STDOUT)

            # Display results after run
            with open(result_file) as f:
                print(f.read())
        else:
            missing = [f for f in [pan_grandparent, pan027, pan028] if not os.path.exists(f)]
            print(f"  Skipping {chrom}: missing {', '.join(missing)}")

    return summary


# ---- Example usage ----
results = process_centromere_transmissions("verkko_2.0_centromere_transmissions_chromosome_transmissions.csv")

print(len(results["two_generational_inheritances"]))
chrom_summary = loop_by_chromosome(results,"/private/groups/migalab/jmmenend/PAN/CDR_realignment_analysis/PAN_centromere_realignment_13101488") #lr_hq
#chrom_summary = loop_by_chromosome(results,"/private/groups/migalab/jmmenend/PAN/CDR_realignment_analysis/PAN_centromere_realignment_13101294") #map_ont
