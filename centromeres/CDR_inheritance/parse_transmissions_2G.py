import os
import glob
import shutil
import subprocess
import pandas as pd


def read_2g_transmissions(txt_path):
    """
    Read assignments from 2G.centromere.transmissions.txt.

    Returns
    -------
    dict
        {
            "chr1": {
                "generation1": "...",
                "generation2": "..."
            },
            ...
        }
    """

    df = pd.read_csv(txt_path, sep="\t")

    chromosomes = list(df.columns[1:])

    assignments = {}

    for chrom in chromosomes:
        assignments[chrom] = {
            "generation1": df.loc[0, chrom],
            "generation2": df.loc[1, chrom]
        }

    return assignments


def loop_by_chromosome(assignments, src_folder, verbose=True):

    summary = {}

    for chrom in assignments:

        generation1 = assignments[chrom]["generation1"]
        generation2 = assignments[chrom]["generation2"]

        summary[chrom] = {
            "generation1": generation1,
            "generation2": generation2
        }

        if verbose:
            print(f"\n=== {chrom} ===")
            print(f"Generation 1 : {generation1}")
            print(f"Generation 2 : {generation2}")

        # ----------------------------------------------------------
        # Copy BED files for this chromosome
        # ----------------------------------------------------------

        copied_beds = []

        patterns = [
            f"{generation1}*{generation2}*centrodip.bed",
            f"{generation2}*{generation2}*centrodip.bed",
        ]

        bed_candidates = []

        for pattern in patterns:
            print(pattern)
            bed_candidates.extend(glob.glob(os.path.join(src_folder, pattern)))

        bed_candidates = sorted(set(bed_candidates))

        print(f"  Bed files: {bed_candidates}")

        if len(bed_candidates) != 2:
            raise ValueError(
                f"Expected exactly 2 BED files for chromosome {chrom}, "
                f"found {len(bed_candidates)}: {bed_candidates}"
            )

        for src in bed_candidates:
            dst = os.path.join(os.getcwd(), os.path.basename(src))
            shutil.copy(src, dst)
            copied_beds.append(dst)

        if not copied_beds:
            print(f"  No BED files found for {chrom}")
            continue

        # ----------------------------------------------------------
        # Convert BED files
        # ----------------------------------------------------------

        for bed in copied_beds:

            base = os.path.basename(bed)

            try:
                sample = base.split(".chr")[0]
                chrom = base.split(".")[1]
                hap = base.split(".realigned_to.")[0].split(".")[-1]
                entity = f"{sample}.{chrom}.{hap}"

            except Exception:
                print(f"Skipping malformed file: {base}")
                continue

            out_bed = f"{entity}.bed"

            cmd = (
                f'egrep --color "{chrom}" "{base}" | '
                f'grep -v "transition_CDR" | '
                f'awk "{{ if ((\\$3 - \\$2) >= 3000) print }}" | '
                f'python keep_main_dbscan_cluster_filter_outliers.py > "{out_bed}"'
            )

            print(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, check=False)

        # ----------------------------------------------------------
        # Compare the two generations
        # ----------------------------------------------------------

        generation1_bed = f"{generation1}.bed"
        generation2_bed = f"{generation2}.bed"

        if all(os.path.exists(f) for f in [generation1_bed, generation2_bed]):

            cmd = (
                'source /opt/miniconda/etc/profile.d/conda.sh && '
                'conda activate /private/groups/migalab/mcechova/conda/genomics-r && '
                f'Rscript compare_gens_overlap_span_2G.R '
                f'"{generation1_bed}" "{generation2_bed}"'
            )

            print(f"Comparing: {cmd}")

            result_file = f"{chrom}_CDR_2G_stats.minL3k.txt"

            with open(result_file, "w") as f:
                subprocess.run(
                    cmd,
                    shell=True,
                    executable="/bin/bash",
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    check=False,
                )

            with open(result_file) as f:
                print(f.read())

        else:
            missing = [
                f
                for f in [generation1_bed, generation2_bed]
                if not os.path.exists(f)
            ]

            print(f"Skipping {chrom}: missing {', '.join(missing)}")

    return summary

assignments = read_2g_transmissions(
    "2G.centromere.transmissions.txt"
)

chrom_summary = loop_by_chromosome(
    assignments,
    "/private/groups/migalab/jmmenend/PAN/CDR_realignment_analysis/PAN_centromere_realignment_13101294"
)

#/private/groups/migalab/jmmenend/PAN/CDR_realignment_analysis/PAN_centromere_realignment_13101488 #lr_hq
#/private/groups/migalab/jmmenend/PAN/CDR_realignment_analysis/PAN_centromere_realignment_13101294 #map_ont