#!/bin/bash

#SBATCH --job-name=pedigree_centromere_realignment
#SBATCH --partition=long
#SBATCH --cpus-per-task=64
#SBATCH --mem=512G
#SBATCH --time=7-00:00:00
#SBATCH --output=pedigree_centromere_realignment.%j.log

THREADS=64

#############################
### === Static Params === ###
#############################

OUTDIR="$(pwd)/PAN_centromere_realignment_${SLURM_JOB_ID}"
mkdir -p "$OUTDIR"
TMPDIR="/data/tmp/jmmenend/PAN_realign"
mkdir -p "$TMPDIR"
TMPFILES=()
trap 'rm -f -- "${TMPFILES[@]}"' EXIT INT TERM HUP

# --- Load in CenSat Annotations
PAN010_HAP1_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN010_hap1_HiFi_element_final_hap1.polished.cenSat.bed"
PAN010_HAP2_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN010_hap2_HiFi_element_final_hap2.polished.cenSat.bed"
PAN010_CENSAT="${TMPDIR}/PAN010_HiFi_element_final.diploid.polished.censat.bed"
cat "$PAN010_HAP1_CENSAT" "$PAN010_HAP2_CENSAT" > "$PAN010_CENSAT"

PAN011_HAP1_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN011_hap1_HiFi_element_final_XY_hap1.polished.cenSat.bed"
PAN011_HAP2_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN011_hap2_HiFi_element_final_XY_hap2.polished.cenSat.bed"
PAN011_CENSAT="${TMPDIR}/PAN011_HiFi_element_final.diploid.polished.censat.bed"
cat "$PAN011_HAP1_CENSAT" "$PAN011_HAP2_CENSAT" > "$PAN011_CENSAT"

PAN027_MAT_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN027_mat_HiFi_element_final_mat.polished.cenSat.bed"
PAN027_PAT_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN027_pat_HiFi_element_final_pat.polished.cenSat.bed"
PAN027_CENSAT="${TMPDIR}/PAN027_HiFi_element_final.diploid.polished.censat.bed"
cat "$PAN027_MAT_CENSAT" "$PAN027_PAT_CENSAT" > "$PAN027_CENSAT"

PAN028_HAP1_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN028_hap1_HiFi_element_final_hap1.polished.cenSat.bed"
PAN028_HAP2_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN028_hap2_HiFi_element_final_hap2.polished.cenSat.bed"
PAN028_CENSAT="${TMPDIR}/PAN028_HiFi_element_final.diploid.polished.censat.bed"
cat "$PAN028_HAP1_CENSAT" "$PAN028_HAP2_CENSAT" > "$PAN028_CENSAT"

BEDS=( 
    "$PAN010_CENSAT" 
    "$PAN011_CENSAT" 
    "$PAN027_CENSAT" 
    "$PAN028_CENSAT" 
)

# --- Load in BAM Files 
PAN010_BAM="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/adaptive_sampling/PAN010.AS.fastq.cpg.winnowmap.assembly.v1.0.PAN010.diploid.bam"
PAN011_BAM="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/adaptive_sampling/PAN011.AS.fastq.cpg.winnowmap.assembly.v1.0.PAN011.diploid.bam"
PAN027_BAM="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/adaptive_sampling/PAN027.AS.fastq.cpg.winnowmap.assembly.v1.0.PAN027.diploid.bam"
PAN028_BAM="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/adaptive_sampling/PAN028.AS.fastq.cpg.winnowmap.assembly.v1.0.PAN028.diploid.bam"
BAMS=( 
    "$PAN010_BAM" 
    "$PAN011_BAM" 
    "$PAN027_BAM" 
    "$PAN028_BAM"
)

# --- Load in FASTA Assemblies 
PAN010_HAP1_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN010_hap1_HiFi_element_final_hap1.polished.fasta"
PAN010_HAP2_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN010_hap2_HiFi_element_final_hap2.polished.fasta"
PAN010_REF="${TMPDIR}/PAN010_HiFi_element_final.diploid.polished.fasta"
cat "$PAN010_HAP1_REF" "$PAN010_HAP2_REF" > "$PAN010_REF"
samtools faidx "$PAN010_REF"

PAN011_HAP1_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN011_hap1_HiFi_element_final_XY_hap1.polished.fasta"
PAN011_HAP2_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN011_hap2_HiFi_element_final_XY_hap2.polished.fasta"
PAN011_REF="${TMPDIR}/PAN011_HiFi_element_final.diploid.polished.fasta"
cat "$PAN011_HAP1_REF" "$PAN011_HAP2_REF" > "$PAN011_REF"
samtools faidx "$PAN011_REF"

PAN027_MAT_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN027_mat_HiFi_element_final_mat.polished.fasta"
PAN027_PAT_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN027_pat_HiFi_element_final_pat.polished.fasta"
PAN027_REF="${TMPDIR}/PAN027_HiFi_element_final.diploid.polished.fasta"
cat "$PAN027_MAT_REF" "$PAN027_PAT_REF" > "$PAN027_REF"
samtools faidx "$PAN027_REF"

PAN028_HAP1_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN028_hap1_HiFi_element_final_hap1.polished.fasta"
PAN028_HAP2_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN028_hap2_HiFi_element_final_hap2.polished.fasta"
PAN028_REF="${TMPDIR}/PAN028_HiFi_element_final.diploid.polished.fasta"
cat "$PAN028_HAP1_REF" "$PAN028_HAP2_REF" > "$PAN028_REF"
samtools faidx "$PAN028_REF"

REFS=( 
    "$PAN010_REF" 
    "$PAN011_REF" 
    "$PAN027_REF" 
    "$PAN028_REF" 
)

SAMPLES=(
    "PAN010"
    "PAN011"
    "PAN027"
    "PAN028"
)

CHROMS=(
    "chr1" "chr2" "chr3" "chr4" "chr5"
    "chr6" "chr7" "chr8" "chr9" "chr10"
    "chr11" "chr12" "chr13" "chr14" "chr15"
    "chr16" "chr17" "chr18" "chr19" "chr20"
    "chr21" "chr22" "chrX" 
)

CENSAT_TERM="active_hor"
MIN_MAPQ=10
MIN_READ_LEN=10000
HAPLOTYPES=("haplotype1" "haplotype2")
PAN027_HAPLOTYPES=("paternal" "maternal")

echo "Setting ${SLURM_JOB_ID} Run Parameters to..."
echo "  MINIMAP2 PRESET = map-ont"
echo "  MIN_MAPQ = $MIN_MAPQ"
echo "  MIN_READ_LEN = $MIN_READ_LEN"
echo "  CENSAT_TERM = $CENSAT_TERM"
echo "  THREADS = $THREADS"

METADATA="pan_centromere_realignment.${SLURM_JOB_ID}.metadata.tsv"
echo -e "from_chrom\tprealign_count\tto_chrom\tpostalign_count\tpercentage" > "$METADATA"

source /private/groups/migalab/jmmenend/.software/anaconda3/etc/profile.d/conda.sh

set -euo pipefail
set -x

###################################
### === Main Loop of Script === ###
###################################
# Loop through each BAM and Censat bed file
#   * Perform alignment to matched assembly
#   Loop through chromosomes
#       Loop through haplotypes 
#           * Extract reads for one chrom/hap combination
#           * Realign those reads to PAN027 (Diploid)
#           Loop through PAN027 Haplotypes of that chromosome
#               * Look at Mapping Efficiency
#               * Call Aggregated Methylation
#               * Call CDRs

conda activate /private/groups/migalab/jmmenend/.conda_envs/minimap2

FQS_FOR_REALIGNMENT=()
for i in "${!BAMS[@]}"; do

    # === Perform Alignment to Matched Diploid Genome (To Fetch Reads) === #

    BAM="${BAMS[$i]}"
    REF="${REFS[$i]}"
    BED="${BEDS[$i]}"
    SAMPLE="${SAMPLES[$i]}"

    TMP_FQ1=$(mktemp /data/tmp/jmmenend/PAN_realign/PAN_realign.${SAMPLE}.XXXXXX.fq); TMPFILES+=("$TMP_FQ1")
    TMP_BAM1=$(mktemp /data/tmp/jmmenend/PAN_realign/PAN_realign.${SAMPLE}.minimap2.XXXXXX..bam); TMPFILES+=("$TMP_BAM1" "${TMP_BAM1}.bai")

    samtools fastq -@ "$THREADS" -T '*' "$BAM" > "$TMP_FQ1"
    minimap2 -ax map-ont -I 16G --eqx --cs -Y -L -y -p 0.5 -k 17 -t "$THREADS" "$REF" "$TMP_FQ1" | \
        samtools view -@ "$THREADS" -q "$MIN_MAPQ" -F 0x904 -bh - | \
        samtools sort -@ "$THREADS" -o "$TMP_BAM1" -
    samtools index -@ "$THREADS" "$TMP_BAM1"

    for CHROM in "${CHROMS[@]}"; do
        for j in "${!HAPLOTYPES[@]}"; do
            if [[ "$SAMPLE" == "PAN027" ]]; then
                HAP="${PAN027_HAPLOTYPES[$j]}"
            else
                HAP="${HAPLOTYPES[$j]}"
            fi

            TMP_BED1=$(mktemp /data/tmp/jmmenend/PAN_realign/PAN_realign.${SAMPLE}.${CHROM}.${HAP}.XXXXXX.bed); TMPFILES+=("$TMP_BED1")
            TMP_FQ2=$(mktemp /data/tmp/jmmenend/PAN_realign/PAN_realign.${SAMPLE}.${CHROM}.${HAP}.XXXXXX.fq); TMPFILES+=("$TMP_FQ2")
            TMP_BAM2="${OUTDIR}/PAN_realign.${SAMPLE}.${CHROM}.${HAP}.bam"

            # pull out the active array for each contig
            awk -v chrom=".${CHROM}." -v hap="$HAP" -v sat="$CENSAT_TERM" \
                'BEGIN { FS = OFS = "\t" } index($1, chrom) && index($1, hap) && index($4, sat) { print }' \
                "$BED" >> "$TMP_BED1"
            if [[ $(wc -l < "$TMP_BED1") -eq 0 ]]; then
                continue
            fi

            # pull out the reads for only one active array
            # pull out reads over 100kb that align well
            samtools view \
                -@ "$THREADS" \
                -bh \
                -q "$MIN_MAPQ" \
                -e "length(seq)>${MIN_READ_LEN}" \
                -F 0x904 \
                -L "$TMP_BED1" \
                "$TMP_BAM1" > "$TMP_BAM2"
            samtools fastq -@ "$THREADS" -T '*' "$TMP_BAM2" > "$TMP_FQ2"

            FROM_CONTIG=$(awk 'NR==1 {print $1}' "$TMP_BED1")
            PREALIGNMENT_READ_COUNT=$(samtools view -c "$TMP_BAM2")

            TMP_BAM3=$(mktemp /data/tmp/jmmenend/PAN_realign/PAN_realign.${SAMPLE}.${CHROM}.${HAP}.to.PAN027.XXXXXX.bam); TMPFILES+=("$TMP_BAM3" "${TMP_BAM3}.bai")

            #######################################################
            # === the realignment portion of the scripting... === #
            #######################################################
            minimap2 -ax map-ont -I 16G --eqx --cs -Y -L -y -p 0.5 -k 17 -t "$THREADS" "$PAN027_REF" "$TMP_FQ2" | \
                samtools view \
                    -@ "$THREADS" \
                    -bh \
                    -q "$MIN_MAPQ" \
                    -e "length(seq)>${MIN_READ_LEN}" \
                    -F 0x904 | \
                        samtools sort -@ "$THREADS" -o "$TMP_BAM3" -
            samtools index -@ "$THREADS" "$TMP_BAM3"

            # Loop through the two haplotypes of PAN027 - to see the alignment efficiency + modkit + CDR's (on each hap individually)
            for PAN027_HAP in "${PAN027_HAPLOTYPES[@]}"; do

                # Fetching CenSat Annotations for an individual Haplotype 
                TMP_PAN027_BED=$(mktemp /data/tmp/jmmenend/PAN_realign/PAN_realign.${SAMPLE}.${CHROM}.${HAP}.to.PAN027.${CHROM}.${PAN027_HAP}.XXXXXX.bed); TMPFILES+=("$TMP_PAN027_BED")
                awk -v chrom=".${CHROM}." -v hap="$PAN027_HAP" -v sat="$CENSAT_TERM" \
                    'BEGIN { FS = OFS = "\t" } index($1, chrom) && index($1, hap) && index($4, sat) { print }' \
                    "$PAN027_CENSAT" >> "$TMP_PAN027_BED"
                if [[ $(wc -l < "$TMP_PAN027_BED") -eq 0 ]]; then
                    continue
                fi

                TMP_BAM4="${OUTDIR}/PAN_realign.${SAMPLE}.${CHROM}.${HAP}.to.PAN027.${CHROM}.${PAN027_HAP}.bam"
                samtools view \
                    -@ "$THREADS" \
                    -bh \
                    -q "$MIN_MAPQ" \
                    -e "length(seq)>${MIN_READ_LEN}" \
                    -F 0x904 \
                    -L "$TMP_PAN027_BED" "$TMP_BAM3" | \
                        samtools sort -@ "$THREADS" -o "$TMP_BAM4" -
                samtools index -@ "$THREADS" "$TMP_BAM4"

                TO_CONTIG=$(awk 'NR==1 {print $1}' "$TMP_PAN027_BED")
                POSTALIGNMENT_READ_COUNT=$(samtools view -c "$TMP_BAM4")

                PERCENT=$(awk -v p="${POSTALIGNMENT_READ_COUNT}" -v q="${PREALIGNMENT_READ_COUNT}" 'BEGIN { if (q > 0) printf "%.2f%%", (p/q)*100; else print "NA" }')
                echo -e "${FROM_CONTIG}\t${PREALIGNMENT_READ_COUNT}\t${TO_CONTIG}\t${POSTALIGNMENT_READ_COUNT}\t${PERCENT}" >> "$METADATA"

                OUTPUT_BASE="${OUTDIR}/${SAMPLE}.${CHROM}.${HAP}.realigned_to.PAN027.${CHROM}.${PAN027_HAP}"

                MODKIT="${OUTPUT_BASE}.modkit.bed"
                conda activate /private/groups/migalab/jmmenend/.conda_envs/modkit
                modkit pileup \
                    -t "$THREADS" \
                    --force-allow-implicit \
                    --cpg \
                    --ref "$PAN027_REF" \
                    --combine-strands \
                    --mod-thresholds m:0.8 --filter-threshold C:0.5 \
                    "$TMP_BAM4" "$MODKIT"
                conda deactivate

                # call CDR's
                CDR="${OUTPUT_BASE}.centrodip.bed"
                conda activate /private/groups/migalab/jmmenend/.conda_envs/centrodip
                centrodip \
                    --label "CDR" \
                    "$MODKIT" \
                    "$TMP_PAN027_BED" \
                    "$CDR"
                conda deactivate
            done
        done
    done
done

conda deactivate

rm -r "$TMPDIR"

echo "Done."