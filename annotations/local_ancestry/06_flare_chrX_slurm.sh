#!/bin/bash

#SBATCH --job-name=flare
#SBATCH -N 1
#SBATCH -A rmccoy22-condo
#SBATCH --time=00:30:00
#SBATCH --partition=shared
#SBATCH --mem=30G
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu
#SBATCH --array=1-4

# Run flare on chrX

#================#
# Set data paths #
#================#

MANIFEST=sample_list.txt
FLARE=/home/syan11/code/flare/flare.jar
KGP=/scratch16/rmccoy22/syan11/washu_local_ancestry/kgp_multiallele_split
WASHU=/scratch16/rmccoy22/syan11/washu_local_ancestry/washu_vcfs
MAP=/scratch16/rmccoy22/syan11/washu_local_ancestry/genetic_map_chm13
EXCLUDE_PATH=/scratch16/rmccoy22/syan11/washu_local_ancestry/exclude_files
OUTDIR=flare

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $MANIFEST`
echo $sample

ml java
ml bcftools

# Get data paths
refvcf=$KGP/1kgp.chm13.chrX.multiallelic_split.vcf.gz
queryvcf=$WASHU/${sample}_chm13_chrX_final.vcf.gz
mapfile=$MAP/chrX.CHM13v2.map


#===============================================#
# chrX non-PAR (generate new demographic model) #
#===============================================#

out=$OUTDIR/${sample}_chm13_chrX_nonpar
# Reference panel of XX samples only
PANEL=/scratch16/rmccoy22/syan11/washu_local_ancestry/afr_eur.XX.panel
# Exclude all variants in PAR1 and PAR2
EXCLUDE=${EXCLUDE_PATH}/$sample.exclude_without_nonpar.txt

java -Xmx25g -jar $FLARE \
	ref=$refvcf \
	ref-panel=$PANEL \
	gt=$queryvcf \
	map=$mapfile \
	out=$out \
	excludemarkers=$EXCLUDE \
	probs=true


#=============================================#
# chrX PARs (use autosomal demographic model) #
#=============================================#

# Use autosomal demographic model
model=$OUTDIR/${sample}_chm13_chr1.model
# Reference panel of all samples
PANEL=/scratch16/rmccoy22/syan11/washu_local_ancestry/afr_eur.panel

### PAR1
out=$OUTDIR/${sample}_chm13_chrX_par1
# Exclude all variants outside of PAR1
EXCLUDE=${EXCLUDE_PATH}/$sample.exclude_without_par1.txt

java -Xmx8g -jar $FLARE \
	ref=$refvcf \
	ref-panel=$PANEL \
	gt=$queryvcf \
	map=$mapfile \
	out=$out \
	excludemarkers=$EXCLUDE \
	probs=true \
	model=$model

### PAR2
out=$OUTDIR/${sample}_chm13_chrX_par2
# Exclude all variants outside of PAR2
EXCLUDE=${EXCLUDE_PATH}/$sample.exclude_without_par2.txt

java -Xmx8g -jar $FLARE \
	ref=$refvcf \
	ref-panel=$PANEL \
	gt=$queryvcf \
	map=$mapfile \
	out=$out \
	excludemarkers=$EXCLUDE \
	probs=true \
	model=$model
