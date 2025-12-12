#!/bin/bash

#SBATCH --job-name=PAN_centrodip
#SBATCH --partition=long
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --output=PAN_centrodip.%j.log

source /private/groups/migalab/jmmenend/.software/anaconda3/etc/profile.d/conda.sh

set -eux -o pipefail

THREADS=32
WORKDIR=$(pwd)

mkdir -p "/data/tmp/jmmenend"
TMPDIR=$(mktemp -d "/data/tmp/jmmenend/PAN027_centrodip_workflow_XXXXXX"); 
cd ${TMPDIR}

trap "rm -rf ${TMPDIR}" EXIT INT TERM HUP

PAN010_HAP1_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN010_hap1_HiFi_element_final_hap1.polished.cenSat.bed"
PAN010_HAP2_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN010_hap2_HiFi_element_final_hap2.polished.cenSat.bed"
PAN010_CENSAT="PAN010_HiFi_element_final.diploid.polished.cenSat.bed"
cat "$PAN010_HAP1_CENSAT" "$PAN010_HAP2_CENSAT" > "$PAN010_CENSAT"

PAN011_HAP1_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN011_hap1_HiFi_element_final_XY_hap1.polished.cenSat.bed"
PAN011_HAP2_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN011_hap2_HiFi_element_final_XY_hap2.polished.cenSat.bed"
PAN011_CENSAT="PAN011_HiFi_element_final.diploid.polished.cenSat.bed"
cat "$PAN011_HAP1_CENSAT" "$PAN011_HAP2_CENSAT" > "$PAN011_CENSAT"

PAN027_MAT_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN027_mat_HiFi_element_final_mat.polished.cenSat.bed"
PAN027_PAT_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN027_pat_HiFi_element_final_pat.polished.cenSat.bed"
PAN027_CENSAT="PAN027_HiFi_element_final.diploid.polished.cenSat.bed"
cat "$PAN027_MAT_CENSAT" "$PAN027_PAT_CENSAT" > "$PAN027_CENSAT"

PAN028_HAP1_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN028_hap1_HiFi_element_final_hap1.polished.cenSat.bed"
PAN028_HAP2_CENSAT="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/PAN028_hap2_HiFi_element_final_hap2.polished.cenSat.bed"
PAN028_CENSAT="PAN028_HiFi_element_final.diploid.polished.cenSat.bed"
cat "$PAN028_HAP1_CENSAT" "$PAN028_HAP2_CENSAT" > "$PAN028_CENSAT"

BEDS=( 
    "$PAN010_CENSAT" 
    "$PAN011_CENSAT" 
    "$PAN027_CENSAT" 
    "$PAN028_CENSAT" 
)

PAN010_BAM="/private/nanopore/basecalled/pedigree/HG06803/06_18_24_R1041_UL_Pedigree_HG06803_1_dorado0.7.2_sup5.0.0_5mCG_5hmCG.bam"
PAN011_BAM="/private/nanopore/basecalled/pedigree/HG06804/06_18_24_R1041_UL_Pedigree_HG06804_1_dorado0.7.2_sup5.0.0_5mCG_5hmCG.bam"
PAN027_BAM="/private/nanopore/basecalled/pedigree/HG06807/06_18_24_R1041_UL_Pedigree_HG06807_1_dorado0.7.2_sup5.0.0_5mCG_5hmCG.bam"
PAN028_BAM="/private/nanopore/basecalled/pedigree/HG06808/06_18_24_R1041_UL_Pedigree_HG06808_1_dorado0.7.2_sup5.0.0_5mCG_5hmCG.bam"
BAMS=( 
    "$PAN010_BAM" 
    "$PAN011_BAM" 
    "$PAN027_BAM" 
    "$PAN028_BAM" 
)

PAN010_HAP1_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN010_hap1_HiFi_element_final_hap1.polished.fasta"
PAN010_HAP2_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN010_hap2_HiFi_element_final_hap2.polished.fasta"
PAN010_REF="PAN010_HiFi_element_final.diploid.polished.fasta"
cat "$PAN010_HAP1_REF" "$PAN010_HAP2_REF" > "$PAN010_REF"
samtools faidx "$PAN010_REF"

PAN011_HAP1_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN011_hap1_HiFi_element_final_XY_hap1.polished.fasta"
PAN011_HAP2_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN011_hap2_HiFi_element_final_XY_hap2.polished.fasta"
PAN011_REF="PAN011_HiFi_element_final.diploid.polished.fasta"
cat "$PAN011_HAP1_REF" "$PAN011_HAP2_REF" > "$PAN011_REF"
samtools faidx "$PAN011_REF"

PAN027_MAT_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN027_mat_HiFi_element_final_mat.polished.fasta"
PAN027_PAT_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN027_pat_HiFi_element_final_pat.polished.fasta"
PAN027_REF="PAN027_HiFi_element_final.diploid.polished.fasta"
cat "$PAN027_MAT_REF" "$PAN027_PAT_REF" > "$PAN027_REF"
samtools faidx "$PAN027_REF"

PAN028_HAP1_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN028_hap1_HiFi_element_final_hap1.polished.fasta"
PAN028_HAP2_REF="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/PAN028_hap2_HiFi_element_final_hap2.polished.fasta"
PAN028_REF="PAN028_HiFi_element_final.diploid.polished.fasta"
cat "$PAN028_HAP1_REF" "$PAN028_HAP2_REF" > "$PAN028_REF"
samtools faidx "$PAN028_REF"

REFS=( 
    "$PAN010_REF" 
    "$PAN011_REF" 
    "$PAN027_REF" 
    "$PAN028_REF" 
)

######################################
### === Main Loop to Call CDRs === ###
######################################

for i in "${!BEDS[@]}"; do
    BED="${BEDS[$i]}"
    BAM="${BAMS[$i]}"
    REF="${REFS[$i]}"

    if [[ ! -f "$BED" ]] || [[ ! -f "$BAM" ]] || [[ ! -f "$REF" ]]; then
        echo "[WARN] Skipping index ${i}: missing BED/BAM/REF" >&2
        continue
    fi

    # convert UBAM to FQ with methylation tags
    conda activate samtools
    TMPFQ="$(mktemp ${TMPDIR}/PAN.centrodip.XXXXXX.fq)"
    samtools fastq -@ "$THREADS" -T '*' "$BAM" > "$TMPFQ"
    conda deactivate

    # realign FQ to polished reference
    echo "[INFO] Realign $(basename "${BAM}") → $(basename "${REF}")"
    conda activate minimap2
    TMPBAM="$(mktemp ${TMPDIR}/PAN.centrodip.XXXXXX.bam)"
    minimap2 -ax map-ont -I 16G --eqx --cs -Y -L -y -p 0.5 -k 17 -t "$THREADS" "$REF" "$TMPFQ" | \
        samtools view -@ "$THREADS" -q 10 -F 0x904 -bh - | \
        samtools sort -@ "$THREADS" -o "$TMPBAM" -
    samtools index -@ "$THREADS" "$TMPBAM"
    conda deactivate

    # call modkit on polished reference
    echo "[INFO] modkit pileup"
    conda activate modkit
    MODKIT="${WORKDIR}/$(basename ${BAM}).pileup.bed"
    modkit pileup  -t "$THREADS" --force-allow-implicit --cpg --ref "$REF" --combine-strands --mod-thresholds m:0.8 --filter-threshold C:0.5 "$TMPBAM" "$MODKIT"
    conda deactivate

    # remove rows with values under 5 in the fourth column
    MODKIT_COV5="${MODKIT%.bed}.cov5.bed"
    awk '$4 >= 5' "$MODKIT" > "$MODKIT_COV5"

    # subset censats to only active alpha
    ACTIVE_HOR="${BED%.bed}.active_hor.bed"
    awk '$4 ~ /active_hor/' "$BED" > "$ACTIVE_HOR"

    # centrodip on modkit file
    echo "[INFO] centrodip"
    conda activate centrodip
    CENTRODIP="${WORKDIR}/$(basename ${BAM}).centrodip.bed"
    centrodip "$MODKIT_COV5" "$ACTIVE_HOR" "$CENTRODIP"
    conda deactivate

done