#!/bin/bash

# ============================================================
# USER CONFIGURATION — edit these paths before running
# ============================================================
PROJECT_DIR="/path/to/your/project"     # root project directory
SCRATCH_DIR="/path/to/your/scratch"     # scratch space for jobstore (~100 GB free needed)
IMAGE="quay.io/comparative-genomics-toolkit/cactus:v3.1.4"  # Cactus Podman image

# Input data directory — FASTAs expected at:
#   DATA_DIR/CHM13_per_chromosome/chr{N}.fasta
#   DATA_DIR/{SAMPLE}/{SAMPLE}.chr{N}.haplotype{1|2}.fasta
#   DATA_DIR/PAN027/PAN027.chr{N}.maternal.fasta
#   DATA_DIR/PAN027/PAN027.chr{N}.paternal.fasta
#   DATA_DIR/PAN028/PAN028.chr{N}.haplotype{1|2}.fa  (note: .fa not .fasta)
DATA_DIR="${PROJECT_DIR}/data"
# ============================================================

RESULTS_DIR="${PROJECT_DIR}/results"
SEQFILE_TMP="${SCRATCH_DIR}/seqfile_tmp.txt"

mkdir -p "$RESULTS_DIR"
mkdir -p "${SCRATCH_DIR}/jobstore"
mkdir -p "${PROJECT_DIR}/logs"

cd "$PROJECT_DIR" || exit 1

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do

    # PAN011_hap2: for chrX use chrY (PAN011 is XY — no chrX hap2 exists)
    if [ "$CHR" = "X" ]; then
        PAN011_HAP2="${DATA_DIR}/PAN011/PAN011.chrY.haplotype2.fasta"
    else
        PAN011_HAP2="${DATA_DIR}/PAN011/PAN011.chr${CHR}.haplotype2.fasta"
    fi

    # Generate per-chromosome seqfile (CHM13 included as a sample)
    cat > "$SEQFILE_TMP" <<EOF
CHM13       ${DATA_DIR}/CHM13_per_chromosome/chr${CHR}.fasta
PAN010_hap1	${DATA_DIR}/PAN010/PAN010.chr${CHR}.haplotype1.fasta
PAN010_hap2	${DATA_DIR}/PAN010/PAN010.chr${CHR}.haplotype2.fasta
PAN011_hap1	${DATA_DIR}/PAN011/PAN011.chr${CHR}.haplotype1.fasta
PAN011_hap2	${PAN011_HAP2}
PAN027_mat	${DATA_DIR}/PAN027/PAN027.chr${CHR}.maternal.fasta
PAN027_pat	${DATA_DIR}/PAN027/PAN027.chr${CHR}.paternal.fasta
PAN028_hap1	${DATA_DIR}/PAN028/PAN028.chr${CHR}.haplotype1.fa
PAN028_hap2	${DATA_DIR}/PAN028/PAN028.chr${CHR}.haplotype2.fa
EOF

    for REF in PAN027_mat PAN027_pat; do

        NAME="PAN027.chr${CHR}.project_${REF}"
        OUT_DIR="${RESULTS_DIR}/${NAME}"
        JOBSTORE="${SCRATCH_DIR}/jobstore/${NAME}"

        mkdir -p "$OUT_DIR"

        if [ -d "$JOBSTORE" ]; then
            echo "Removing old jobstore for ${NAME}..."
            rm -rf "$JOBSTORE"
        fi

        echo "===== Running ${NAME} ====="

        nice -n 15 podman run --rm \
            --cgroup-manager=cgroupfs \
            --userns=keep-id:uid=1000,gid=1000 \
            -v "${PROJECT_DIR}:${PROJECT_DIR}:Z" \
            -v "${SCRATCH_DIR}:${SCRATCH_DIR}:Z" \
            -w "${PROJECT_DIR}" \
            "$IMAGE" \
            cactus-pangenome "${JOBSTORE}" "$SEQFILE_TMP" \
            --outName "$NAME" \
            --outDir "$OUT_DIR" \
            --reference "$REF" \
            --vcf --vcfwave --gbz --chrom-vg \
            --batchSystem single_machine \
            --maxCores 64 \
            --binariesMode local \
            --logFile "logs/${NAME}.log"

        if [ $? -eq 0 ]; then
            echo "${NAME} done. Removing jobstore..."
            rm -rf "$JOBSTORE"
        else
            echo "ERROR in ${NAME}. Jobstore kept at ${JOBSTORE} for inspection."
        fi

    done

done

echo "===== All chromosomes done ====="
