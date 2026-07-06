#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --job-name="align_centromeres.20250527"
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --partition=medium
#SBATCH --output=align_centromeres.20250527.%j.log

set -e
set -x

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/.conda/envs/methylation

numCPU=64

input_fasta="$1"

# Check for input
if [ -z "$1" ]; then
  echo "❌ Error: No input FASTA file provided."
  echo "Usage: align_centromeres <input.fa>"
  exit 1
fi

base_name=$(basename "$input_fasta")

awk -v prefix="$base_name" '
  /^>/ {n++}
  {
    if (n==1) out=prefix "_grandparent.fa";
    else if (n==2) out=prefix "_mother.fa";
    else if (n==3) out=prefix "_granddaughter.fa";
    print > out
  }
' "$input_fasta"

# Align the sequences using minimap2

output_dir="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/after_polishing_20250524/outputs/fasta"

alignment_file="$output_dir/${base_name}_generation1_aligned.paf"
stat_file="$output_dir/${base_name}_generation1_stat.txt"
vcf_file_original="$output_dir/${base_name}_generation1_original.vcf"
vcf_file_shifted="$output_dir/${base_name}_generations1_shifted.vcf"

minimap2 -t ${numCPU} --cs -x asm5 -c "${input_fasta}_mother.fa" "${input_fasta}_grandparent.fa" > "$alignment_file" 2> /dev/null
paftools.js stat "$alignment_file" >${stat_file}
paftools.js call -f "${input_fasta}_mother.fa" "$alignment_file" > "$vcf_file_original"


#Preserve header lines (those starting with #), Process only data lines, Output a valid, shifted VCF
awk -F'\t' 'BEGIN { OFS = "\t" }
    /^#/ { print; next }  # keep header lines unmodified
    {
        split($1, coords, /[:\-]/);
        $2 = $2 + coords[2] - 1;
        $1 = coords[1];
        print
    }
' "${vcf_file_original}" > "${vcf_file_shifted}"

#rm $alignment_file

alignment_file="$output_dir/${input_fasta}_generation2_aligned.paf"
stat_file="$output_dir/${input_fasta}_generation2_stat.txt"
vcf_file_original="$output_dir/${input_fasta}_generation2_original.vcf"
vcf_file_shifted="$output_dir/${base_name}_generations2_shifted.vcf"

minimap2 -t ${numCPU} --cs -x asm5 -c "${input_fasta}_mother.fa" "${input_fasta}_granddaughter.fa" > "$alignment_file" 2> /dev/null
paftools.js stat "$alignment_file" >${stat_file}
paftools.js call -f "${input_fasta}_mother.fa" "$alignment_file" > "$vcf_file_original"

#Preserve header lines (those starting with #), Process only data lines, Output a valid, shifted VCF
awk -F'\t' 'BEGIN { OFS = "\t" }
    /^#/ { print; next }  # keep header lines unmodified
    {
        split($1, coords, /[:\-]/);
        $2 = $2 + coords[2] - 1;
        $1 = coords[1];
        print
    }
' "${vcf_file_original}" > "${vcf_file_shifted}"

#rm $alignment_file

echo "Done."



