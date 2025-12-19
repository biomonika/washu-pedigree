#!/bin/bash

# Make file of excluded sites for flare, removing:
# - Sites not in the dipcall dip.bed files
# - PAR or non-PAR sites for chrX

#=====================#
# Set input variables #
#=====================#

vcfDir=/scratch16/rmccoy22/syan11/washu_local_ancestry/washu_vcfs
kgpDir=/scratch16/rmccoy22/syan11/washu_local_ancestry/kgp_nogenotypes
bedDir=/scratch16/rmccoy22/syan11/washu_local_ancestry/raw_data
excludeDir=/scratch16/rmccoy22/syan11/washu_local_ancestry/exclude_files
sample=$1

ml bedtools
ml bcftools


#===================================================#
# Remove sites not in dipcall confident bed regions #
#===================================================#

bed=$bedDir/assembly.v1.0.$sample.dipcall.bed
tmp=$sample.exclude.tmp.txt
out=$excludeDir/$sample.exclude.txt

for i in {1..22} X
do
	vcfWashu=$vcfDir/${sample}_chm13_chr${i}_final.vcf.gz
	bedtools subtract \
		-a $vcfWashu \
		-b $bed | \
		cut -f 1-2 | grep -v "#" | tr "\t" ":" \
		>> $tmp

	vcfKgp=$kgpDir/1kgp_chm13_chr${i}_nogenotypes.vcf.gz
	bcftools view \
		-i 'N_ALT>1' \
		$vcfKgp | \
		cut -f 1-2 | grep -v "#" | tr "\t" ":" \
		>> $tmp
done

# Sort variants by chromosome position
sort -V $tmp > $out
rm $tmp


#====================================#
# Make exclude files for chrX nonPAR #
#====================================#

kgpX=$kgpDir/1kgp_chm13_chrX_nogenotypes.vcf.gz
nonPar=$excludeDir/chrX_without_nonpar.txt

# Select variants to exclude in non-PAR region
bcftools view \
    -t ^chrX:2394410-153925834 \
    $kgpX | \
    cut -f 1-2 | grep -v "#" | tr "\t" ":" \
    > $nonPar

# Add variants from other regions to exclude file
cat $nonPar $out | \
	sort -V \
	> ${excludeDir}/$sample.exclude_without_nonpar.txt


#==================================#
# Make exclude files for chrX PARs #
#==================================#

kgpX=$kgpDir/1kgp_chm13_chrX_nogenotypes.vcf.gz
par1=$excludeDir/chrX_without_par1.txt
par2=$excludeDir/chrX_without_par2.txt

# Select variants to exclude for PAR1 region
bcftools view \
    -t ^chrX:1-2394410 \
    $kgpX | \
    cut -f 1-2 | grep -v "#" | tr "\t" ":" \
    > $par1
# Select variants to exclude for PAR2 region
bcftools view \
    -t chrX:1-153925834 \
    $kgpX | \
    cut -f 1-2 | grep -v "#" | tr "\t" ":" \
    > $par2

# Add variants from other regions to exclude files
cat $par1 $out | \
	sort -V \
	> ${excludeDir}/$sample.exclude_without_par1.txt
cat $par2 $out | \
	sort -V \
	> ${excludeDir}/$sample.exclude_without_par2.txt