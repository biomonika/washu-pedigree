# Shared blocks represent blocks of DNA transmitted through all three generations
# Converting v1.0 coordinates into v1.1 coordinates and filtering (must contain matching chromosomes and the correct maternal/paternal lineage)

##############################################
# convert from shared blocks to bed files so that we can project/liftover
##############################################

sbatch project_shared_contigs_to_bed_files.sh shared_contigs_maternal_coordinates.tsv
sbatch project_shared_contigs_to_bed_files.sh shared_contigs_paternal_coordinates.tsv

##############################################
# split blocks by haplotype
##############################################

cat maternal.PAN010.blocks.bed | grep "haplotype1" | sort -gk4 >maternal.PAN010.blocks.haplotype1.bed
cat maternal.PAN010.blocks.bed | grep "haplotype2" | sort -gk4 >maternal.PAN010.blocks.haplotype2.bed
cat maternal.PAN027.blocks.bed | grep "maternal" | sort -gk4 >maternal.PAN027.blocks.maternal.bed
cat maternal.PAN028.blocks.bed | grep "haplotype1" | sort -gk4 >maternal.PAN028.blocks.haplotype1.bed
cat maternal.PAN028.blocks.bed | grep "haplotype2" | sort -gk4 >maternal.PAN028.blocks.haplotype2.bed

cat paternal.PAN011.blocks.bed | grep "haplotype1" | sort -gk4 >paternal.PAN011.blocks.haplotype1.bed
cat paternal.PAN011.blocks.bed | grep "haplotype2" | sort -gk4 >paternal.PAN011.blocks.haplotype2.bed
cat paternal.PAN027.blocks.bed | grep "paternal" | sort -gk4 >paternal.PAN027.blocks.paternal.bed
cat paternal.PAN028.blocks.bed | grep "haplotype1" | sort -gk4 >paternal.PAN028.blocks.haplotype1.bed
cat paternal.PAN028.blocks.bed | grep "haplotype2" | sort -gk4 >paternal.PAN028.blocks.haplotype2.bed

##############################################
#liftover the coordinates
##############################################

#PAN010
sbatch lift_coordinates.sh assembly.v1.0.PAN010.diploid_haplotype1.fasta PAN010_hap1_HiFi_element_final_hap1.polished.fasta maternal.PAN010.blocks.haplotype1.bed
sbatch lift_coordinates.sh assembly.v1.0.PAN010.diploid_haplotype2.fasta PAN010_hap2_HiFi_element_final_hap2.polished.fasta maternal.PAN010.blocks.haplotype2.bed

#PAN011
sbatch lift_coordinates.sh assembly.v1.0.PAN011.diploid_haplotype1.fasta PAN011_hap1_HiFi_element_final_hap1.polished.fasta paternal.PAN011.blocks.haplotype1.bed
sbatch lift_coordinates.sh assembly.v1.0.PAN011.diploid_haplotype2.fasta PAN011_hap2_HiFi_element_final_hap2.polished.fasta paternal.PAN011.blocks.haplotype2.bed

#PAN027
sbatch lift_coordinates.sh assembly.v1.0.PAN027.diploid_maternal.fasta PAN027_mat_HiFi_element_final_mat.polished.fasta maternal.PAN027.blocks.maternal.bed
sbatch lift_coordinates.sh assembly.v1.0.PAN027.diploid_paternal.fasta PAN027_pat_HiFi_element_final_pat.polished.fasta paternal.PAN027.blocks.paternal.bed

#PAN028
sbatch lift_coordinates.sh assembly.v1.0.PAN028.diploid_haplotype1.fasta PAN028_hap1_HiFi_element_final_hap1.polished.fasta maternal.PAN028.blocks.haplotype1.bed
sbatch lift_coordinates.sh assembly.v1.0.PAN028.diploid_haplotype2.fasta PAN028_hap2_HiFi_element_final_hap2.polished.fasta maternal.PAN028.blocks.haplotype2.bed

sbatch lift_coordinates.sh assembly.v1.0.PAN028.diploid_haplotype1.fasta PAN028_hap1_HiFi_element_final_hap1.polished.fasta paternal.PAN028.blocks.haplotype1.bed
sbatch lift_coordinates.sh assembly.v1.0.PAN028.diploid_haplotype2.fasta PAN028_hap2_HiFi_element_final_hap2.polished.fasta paternal.PAN028.blocks.haplotype2.bed


##############################################
#combine lifted coordinates back to the file with shared blocks
##############################################

python combine_liftover_back_to_maternal_blocks.py 
python combine_liftover_back_to_paternal_blocks.py 