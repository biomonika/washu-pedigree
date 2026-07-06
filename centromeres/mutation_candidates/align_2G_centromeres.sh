#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --job-name="align_2G_centromeres.20260606"
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --partition=medium
#SBATCH --output=align_centromeres.20260606.%j.log

set -e
set -x

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/.conda/envs/methylation

numCPU=64

chromosome="$1"
base_name=$(basename "$chromosome")
echo "Chromosome: " + $chromosome

# Align the sequences using minimap2
output_dir="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/2G/alignments"

alignment_file="$output_dir/${base_name}_generation1_aligned.paf"
stat_file="$output_dir/${base_name}_generation1_stat.txt"
vcf_file_original="$output_dir/${base_name}_generation1_original.vcf"
vcf_file_shifted="$output_dir/${base_name}_generations1_shifted.vcf"

minimap2 -t ${numCPU} --cs -x asm5 -c "${chromosome}_mother.fa" "${chromosome}_grandparent.fa" > "$alignment_file" 2> /dev/null
paftools.js stat "$alignment_file" >${stat_file}
paftools.js call -f "${chromosome}_mother.fa" "$alignment_file" > "$vcf_file_original"


#Preserve header lines (those starting with #), Process only data lines, Output a valid, shifted VCF
awk -F'\t' 'BEGIN { OFS = "\t" }
    /^#/ { print; next }  # keep header lines unmodified
    {
        sub(/^::/, "", $1)   # remove leading "::" from first field
        split($1, coords, /[:\-]/);
        $2 = $2 + coords[2] - 1;
        $1 = coords[1];
        print
    }
' "${vcf_file_original}" > "${vcf_file_shifted}"

#rm $alignment_file

echo "Done."



