#!/bin/bash

#SBATCH --job-name=vcf_merge
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --partition=shared
#SBATCH --mem=5G
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu
#SBATCH --array=1-23

#============#
# Data setup #
#============#

ml bcftools
ml htslib

# Name of sample to run on
sample=$1

# Chromosome name
i=${SLURM_ARRAY_TASK_ID}
if [[ $SLURM_ARRAY_TASK_ID -eq 23 ]]; then
	i=X
fi

# Directory of 1KGP VCF sites with samples removed
SUBSET=kgp_nogenotypes
# Directory of raw dipcall VCFs
RAW_WASHU=/scratch16/rmccoy22/syan11/washu_local_ancestry/raw_data
# Output directory for processed dipcall VCFs
WASHU=/scratch16/rmccoy22/syan11/washu_local_ancestry/washu_vcfs
# Script for editing PAN011 (XY) chrX genotypes into two identical X haplotypes
editPAN011Script=/scratch16/rmccoy22/syan11/washu_local_ancestry/edit_PAN011_vcf.R


#============================#
# 1. Merge WashU and KGP VCF #
#============================#

kgpNoGT=$SUBSET/1kgp_chm13_chr${i}_nogenotypes.vcf.gz
dipcall=${RAW_WASHU}/assembly.v1.0.${sample}.dipcall.vcf.gz
washuMerge=tmp_${sample}_chr${i}_merge.vcf.gz

# Set missing 1KGP sites to ref in WashU sample
# Do not create multiallelic sites (output as separate sites)
if [[ "$sample" == "PAN011" && "$i" == "X" ]]; then
	tmp0=tmp0_$washuMerge
	bcftools merge \
		--missing-to-ref \
		-O z \
		-o $tmp0 \
		-m none \
		$kgpNoGT \
		$dipcall
	# For PAN011 (XY), run a separate script to set half-missing chrX genotypes to homozygous
	$editPAN011Script \
		--input $tmp0 \
		--output $washuMerge
	# bgzip vcf
	bgzip $washuMerge
else
	bcftools merge \
	-O z \
	-o $washuMerge \
	--missing-to-ref \
	-m none \
	$kgpNoGT \
	$dipcall
fi
tabix $washuMerge


#=======================================#
# 2. Intersect to just SNPs in 1KGP VCF #
#=======================================#

washuIsec=tmp_${sample}_chr${i}_isec.vcf.gz

bcftools isec \
	-O z \
	-o $washuIsec \
	-n =2 \
	-w 1 \
	-r chr${i} \
	$washuMerge \
	$kgpNoGT
tabix $washuIsec


#======================#
# 3. Genotype cleanups #
#======================#

washuFinal=$WASHU/${sample}_chm13_chr${i}_final.vcf.gz

# 4. Keep dipcall sample only
# 5. Remove variants with missing genotypes and multi-allelic sites
# 6. Phase `--missing-to-ref` 0/0 genotypes
bcftools view \
	-e 'GT="mis"' \
	-s syndip \
	$washuIsec | \
	bcftools +/home/syan11/code/bcftools-1.18/plugins/setGT.so \
		-O z \
		-o $washuFinal \
		-- \
		-t a \
		-n p
tabix $washuFinal

# remove tmp files
rm tmp_${sample}_chr${i}_*