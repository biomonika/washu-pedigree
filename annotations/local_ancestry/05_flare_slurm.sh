#!/bin/bash

#SBATCH --job-name=flare
#SBATCH -N 1
#SBATCH --time=01:30:00
#SBATCH --partition=shared
#SBATCH --mem=30G
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu
#SBATCH --array=1-4

# Run flare on autosomes

#================#
# Set data paths #
#================#

ml java

MANIFEST=sample_list.txt
FLARE=/home/syan11/code/flare/flare.jar
KGP=/scratch16/rmccoy22/syan11/washu_local_ancestry/kgp_multiallele_split
PANEL=/scratch16/rmccoy22/syan11/washu_local_ancestry/afr_eur.panel
WASHU=/scratch16/rmccoy22/syan11/washu_local_ancestry/washu_vcfs
MAP=/scratch16/rmccoy22/syan11/washu_local_ancestry/genetic_map_chm13
EXCLUDE_PATH=/scratch16/rmccoy22/syan11/washu_local_ancestry/exclude_files
OUTDIR=flare

# Get sample ID from manifest
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $MANIFEST`
exclude=${EXCLUDE_PATH}/$sample.exclude.txt
echo $sample


#================================================#
# Run flare on chr1 to generate population model #
#================================================#

# Get data paths
refvcf=$KGP/1kgp.chm13.chr1.multiallelic_split.vcf.gz
queryvcf=$WASHU/${sample}_chm13_chr1_final.vcf.gz
mapfile=$MAP/chr1.CHM13v2.map
out=$OUTDIR/${sample}_chm13_chr1

# Run flare
java -Xmx25g -jar $FLARE \
	ref=$refvcf \
	ref-panel=$PANEL \
	gt=$queryvcf \
	map=$mapfile \
	out=$out \
	excludemarkers=$exclude \
	probs=true


#==================================================#
# Use chr1 model to run flare on other chromosomes #
#==================================================#

model=$OUTDIR/${sample}_chm13_chr1.model

for i in {2..22} X
do
	refvcf=$KGP/1kgp.chm13.chr$i.multiallelic_split.vcf.gz
	queryvcf=$WASHU/${sample}_chm13_chr${i}_final.vcf.gz
	mapfile=$MAP/chr$i.CHM13v2.map
	out=$OUTDIR/${sample}_chm13_chr$i

	# For all samples
	java -Xmx8g -jar $FLARE \
		ref=$refvcf \
		ref-panel=$PANEL \
		gt=$queryvcf \
		map=$mapfile \
		out=$out \
		excludemarkers=$exclude \
		probs=true \
		model=$model

	# For PAN010, we did not reuse the chr1 model (fails to converge for other chrs)
	# java -Xmx8g -jar $FLARE \
	# 	ref=$refvcf \
	# 	ref-panel=$PANEL \
	# 	gt=$queryvcf \
	# 	map=$mapfile \
	# 	out=$out \
	# 	excludemarkers=$exclude \
	# 	probs=true
done