#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --job-name="lift_coordinates.20250303"
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=48
#SBATCH --time=18:00:00
#SBATCH --mem=64G
#SBATCH --partition=long
#SBATCH --output=lift_coordinates.20250303.%j.log

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/microsatellite

# Check if exactly 3 arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Error: Expected 3 arguments, but got $#."
    echo "Usage: $0 <original_assembly> <polished_assembly> <original_bed_file>"
    exit 1
fi

# Assign arguments to variables
original_assembly=$1
polished_assembly=$2
original_bed_file=$3

# Check if any of the input files are empty
if [ ! -s "$original_assembly" ]; then
    echo "Error: The file '$original_assembly' is empty or does not exist."
    exit 1
fi
if [ ! -s "$polished_assembly" ]; then
    echo "Error: The file '$polished_assembly' is empty or does not exist."
    exit 1
fi
if [ ! -s "$original_bed_file" ]; then
    echo "Error: The file '$original_bed_file' is empty or does not exist."
    exit 1
fi

# If all checks pass, continue with the script execution
echo "All input files are valid. Proceeding with the script..."

# Remove .fa or .fasta extension and extract basenames
original_basename=$(basename "$original_assembly" .fasta)
original_basename=$(basename "$original_basename" .fa)

polished_basename=$(basename "$polished_assembly" .fasta)
polished_basename=$(basename "$polished_basename" .fa)

alignment_file="${original_basename}_${polished_basename}.paf"
filtered_alignment_file="${original_basename}_${polished_basename}.filtered.paf"
output_file="${original_bed_file}.lifted"

numCPU=48

# Check if the alignment exists
if [ -f "$filtered_alignment_file" ]; then
    echo "File $filtered_alignment_file exists. Skipping minimap2 part of the script."
else
    echo "File does not exist. Running the part of the script that depends on the file."
    #minimap2 -ax asm5 ref.fa asm.fa > aln.sam       # assembly to assembly/ref alignment
    minimap2 -t ${numCPU} --cs -x asm5 -c "$original_assembly" "$polished_assembly" > "$alignment_file" 2> /dev/null
    cat "$alignment_file" | awk '{if ($1==$6) print;}' >"$filtered_alignment_file" #only keep matching chromosomes
fi

working_dir="/private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/polishing/polished_by_Mira"

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/groups/migalab/mcechova/conda/convertor

paf2chain -i ${working_dir}/${filtered_alignment_file} >${working_dir}/${filtered_alignment_file}.chain
liftOver ${working_dir}/${original_bed_file} ${working_dir}/${filtered_alignment_file}.chain ${working_dir}/${output_file}.lifted.bed ${working_dir}/${output_file}.unMapped.bed

cat ${working_dir}/${output_file}.lifted.bed | grep -v "#" | wc -l
cat ${working_dir}/${output_file}.unMapped.bed | grep -v "#" | wc -l

echo "Projection completed. Done."
