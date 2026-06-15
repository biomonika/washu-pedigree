import csv

def process_file(file_name):
    """
    Reads a tab-separated BED file and returns a dictionary {id: concatenated string}.
    """
    data = {}
    with open(file_name, mode='r', newline='\n') as file:
        reader = csv.reader(file, delimiter='\t')  # Use tab as delimiter
        for row in reader:
            if len(row) < 4:
                continue  # Skip malformed lines
            concat_str = f"{row[0]}:{row[1]}-{row[2]}"
            data[row[3]] = concat_str  # row[3] is the ID in the fourth column
    return data

def match_and_create_output(pan027_file, PAN011_haplotype1_file, PAN011_haplotype2_file, PAN028_haplotype1_file, PAN028_haplotype2_file, output_file):
    # Load data from the PAN027 file
    pan027_data = process_file(pan027_file)
    first_dict = process_file(PAN011_haplotype1_file)
    second_dict = process_file(PAN011_haplotype2_file)
    third_dict = process_file(PAN028_haplotype1_file)
    fourth_dict = process_file(PAN028_haplotype2_file)

    # Open the output file for writing
    with open(output_file, mode='w', newline='\n') as output:
        writer = csv.writer(output, delimiter='\t', lineterminator='\n')  # Write tab-separated output

        # Iterate through the PAN027 data and create the output rows
        for pan027_id, middle_col in pan027_data.items():
            # Select first_col from first_dict or second_dict based on available ID
            first_col = first_dict.get(pan027_id) or second_dict.get(pan027_id, "")
            # Select third_col from third_dict or fourth_dict based on available ID
            third_col = third_dict.get(pan027_id) or fourth_dict.get(pan027_id, "")

            # Write a row only if middle_col and third_col are not empty
            if middle_col and third_col:
                writer.writerow([first_col, middle_col, third_col])

if __name__ == "__main__":
    # Input files (adjust file names as needed)
    pan027_file = "paternal.PAN027.blocks.paternal.bed.lifted.lifted.bed"
    PAN011_haplotype1_file = "paternal.PAN011.blocks.haplotype1.bed.lifted.lifted.bed"
    PAN011_haplotype2_file = "paternal.PAN011.blocks.haplotype2.bed.lifted.lifted.bed"
    PAN028_haplotype1_file = "paternal.PAN028.blocks.haplotype1.bed.lifted.lifted.bed"
    PAN028_haplotype2_file = "paternal.PAN028.blocks.haplotype2.bed.lifted.lifted.bed"

    # Output file
    output_file = "paternal.blocks.projected.tsv"

    # Run the matching and output creation
    match_and_create_output(pan027_file, PAN011_haplotype1_file, PAN011_haplotype2_file, PAN028_haplotype1_file, PAN028_haplotype2_file, output_file)

    print("Processing complete! Output written to", output_file)
