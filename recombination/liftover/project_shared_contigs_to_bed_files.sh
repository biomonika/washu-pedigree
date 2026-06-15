#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --job-name="project_shared_contigs_to_bed_files.20250303"
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=32
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --partition=long
#SBATCH --output=project_shared_contigs_to_bed_files.20250303.%j.log

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/conda/microsatellite

input_file=$1
row_number=0


declare -A output_files

echo "Processing input file: $input_file"

# Extract lineage from filename (maternal or paternal)
if [[ "$input_file" =~ (maternal|paternal) ]]; then
    lineage=${BASH_REMATCH[1]}
else
    lineage="unknown"
fi

while IFS=$'\t' read -r col1 col2 col3 col4; do
    ((row_number++))
    
    # Ensure the row has exactly 4 columns
    if [ -z "$col4" ]; then
        continue
    fi
    
    # Function to process a column into BED format
    process_column() {
        local col=$1
        local row_num=$2
        
        prefix=$(echo "$col" | cut -d'.' -f1)
        prefix="${lineage}.${prefix}.blocks"

        chr=$(echo "$col" | cut -d':' -f1)
        start=$(echo "$col" | cut -d':' -f2 | cut -d'-' -f1)
        end=$(echo "$col" | cut -d':' -f2 | cut -d'-' -f2)
        
        # Define output file based on prefix before first dot
        if [ -z "${output_files[$prefix]}" ]; then
            output_files[$prefix]="${prefix}.bed"
            echo -n > "${output_files[$prefix]}" # Clear previous file if exists
        fi
        
        echo -e "$chr\t$start\t$end\t$row_num" >> "${output_files[$prefix]}"
    }
    
    # Process and append to respective BED files
    process_column "$col2" "$row_number"
    process_column "$col3" "$row_number"
    process_column "$col4" "$row_number"

done < "$input_file"

echo "Processing complete. BED files created: ${!output_files[@]}"

echo "Done."