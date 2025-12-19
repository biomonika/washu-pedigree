#!/bin/bash

#SBATCH --job-name=multi_split
#SBATCH -N 1
#SBATCH --time=3:00:00
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

# Chromosome name
i=${SLURM_ARRAY_TASK_ID}
if [[ $SLURM_ARRAY_TASK_ID -eq 23 ]]; then
	i=X
fi

# Directory of 1000 Genomes phased VCFs on CHM13v2.0
KGP=/scratch4/rmccoy22/sharedData/populationDatasets/1KGP_NYGC/CHM13v2.0_phased_vcfs
# Output directory for multiallelic split VCFs
OUTDIR=kgp_multiallele_split
# Output directory for VCF of 1KGP sites (all samples removed, for merging with dipcall samples)
SUBSET=kgp_nogenotypes


#==============================================================#
# 1. Split multiallelic and removed unphased sites in 1KGP VCF #
#==============================================================#

out=$OUTDIR/1kgp.chm13.chr$i.multiallelic_split.vcf.gz
kgpVcf=$KGP/1KGP.CHM13v2.0.chr$i.recalibrated.snp_indel.pass.phased.vcf.gz

bcftools norm \
	-m-both \
	$kgpVcf | \
	bcftools view \
		-O z \
		-o $out \
		--phased
tabix $out


#==================================================#
# 2. Make VCF with all samples besides one removed #
#==================================================#

kgpNoGT=$SUBSET/1kgp_chm13_chr${i}_nogenotypes.vcf.gz

bcftools view \
	-O z \
	-o $kgpNoGT \
	-s HG01880 \
	$out
tabix $kgpNoGT
