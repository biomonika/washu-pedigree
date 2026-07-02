#!/bin/bash

# ============================================================
# USER CONFIGURATION — edit these paths before running
# ============================================================
PROJECT_DIR="/path/to/project"                              # root project directory
TEMPLATE="${PROJECT_DIR}/analyse_block.paternal.sh"
OUT_BASE="${PROJECT_DIR}/results_final_filter_mother_gp"    # output folder for filtered VCFs
STAGING_DIR="${PROJECT_DIR}/results_filtered.based_on_blocks.PAN027only"  # temporary staging folder
# ============================================================

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do

    echo "============================================================"
    echo "Running paternal chr${CHR}..."
    echo "============================================================"

    bash "$TEMPLATE" "chr${CHR}"
    EXIT_CODE=$?

    if [ $EXIT_CODE -ne 0 ]; then
        echo "ERROR on chr${CHR} (exit code $EXIT_CODE). Continuing..."
    else
        echo "chr${CHR} done. Moving VCF files..."
        CHR_DIR="${OUT_BASE}/chr${CHR}_pat"
        mkdir -p "$CHR_DIR"
        mv "${STAGING_DIR}/PAN027.chr${CHR}.project_PAN027_pat.wave.filtered.vcf"   "$CHR_DIR/" 2>/dev/null
        mv "${STAGING_DIR}/PAN027.chr${CHR}.project_PAN027_pat.wave.unfiltered.vcf" "$CHR_DIR/" 2>/dev/null
    fi

done

echo "============================================================"
echo "All chromosomes processed."
echo "============================================================"
