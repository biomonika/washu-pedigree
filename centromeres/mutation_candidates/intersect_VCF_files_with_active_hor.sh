#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --job-name="intersect_VCF_files_with_active_hor.20250529"
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --partition=medium
#SBATCH --output=intersect_VCF_files_with_active_hor.20250529.%j.log

#Summary of what this script does:
#Loops over both maternal and paternal
#Picks the corresponding BED and transmissions file.
#For each chromosome, checks both generation VCFs.
#Intersects with the BED file using bedtools.

set -e
set -x

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/.conda/envs/methylation

# Define transmission sources and their respective BED files
declare -A bedfiles=(
    [maternal]="../PAN027_mat_HiFi_element_final_mat.polished.cenSat.active_hor.bed"
    [paternal]="../PAN027_pat_HiFi_element_final_pat.polished.cenSat.active_hor.bed"
)

for parent in maternal paternal; do
    bedfile="${bedfiles[$parent]}"
    transmissions_file="${parent}.transmissions.txt"

    while read chr; do
        for gen in generations1 generations2; do
            vcf_file="${chr}_extracted.fasta_${gen}_shifted.vcf"

            if [[ -f "$vcf_file" ]]; then
                vcf_base="${vcf_file%.*}"
                output="${vcf_base}_${parent}_active_hor.vcf"
                bedtools intersect -header -a "$vcf_file" -b "$bedfile" > "$output"
                echo "Created: $output"
            else
                echo "Missing: $vcf_file" >&2
            fi
        done
    done < "$transmissions_file"
done
