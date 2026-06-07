#!/bin/bash

# ============================================================
# USER CONFIGURATION — edit these paths before running
# ============================================================
SCRIPTS_DIR="/path/to/final_scripts"             # folder where make_stats.per_chr.py and create_chrom_map.py live
OUT_BASE="/path/to/results_filtered_all"         # results folder to process
# ============================================================

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do
    for HAP in mat pat; do
        STATS_DIR="${OUT_BASE}/chr${CHR}_${HAP}"
        if [ ! -d "$STATS_DIR" ]; then
            echo "Skipping chr${CHR}_${HAP} (folder does not exist)"
            continue
        fi
        echo "===== chr${CHR}_${HAP} ====="
        python3 "${SCRIPTS_DIR}/make_stats.per_chr.py" "$CHR" "$HAP" "$OUT_BASE"
        python3 "${SCRIPTS_DIR}/create_chrom_map.py"   "$CHR" "$HAP" "$OUT_BASE"
    done
done

echo "===== Done ====="
