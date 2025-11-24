import sys

def check_sequences(sequence_name, seq1_data, seq2_data, found_start, buffer):
    """
    Compare two aligned sequence chunks and collect mismatch/gap information.

    Parameters
    ----------
    sequence_name : str
        Name of the sequence (used only in output).
    seq1_data, seq2_data : [str, int]
        [0] sequence chunk for sequence 1/2
        [1] current position in sequence 1/2
    found_start : bool
        Whether alignment has yet encountered first A/C/G/T matched base pair.
    buffer : list
        Temporary list to hold mismatch records until validated.

    Returns
    -------
    output : str
        Lines to write to output file.
    line_len : int
        Number of validated (non-gap) positions processed.
    line_mismatches : int
        Number of mismatches in this chunk.
    line_gaps : int
        Number of gap-related mismatches.
    found_start : bool
        Updated state of found_start.
    buffer : list
        Possibly cleared mismatch buffer.
    p1, p2 : int
        Updated positions of seq1 and seq2.
    """

    # seq<N>_data[1] = start position, seq<N>_data[2] = sequence chunk, seq<N>_data[3] = end position
    output = ""
    p1, p2 = seq1_data[1], seq2_data[1]
    seq1_chunk, seq2_chunk = seq1_data[0], seq2_data[0]

    line_len, line_len_curr = 0, 0
    line_mismatches, line_gaps = 0, 0
    for a, b in zip(seq1_chunk, seq2_chunk):
        if not found_start:
            if a in ("A", "C", "G", "T") and b in ("A", "C", "G", "T"):
                found_start = True
            else:
                continue

        line_len_curr += 1

        if a != "-":
            p1 += 1
        if b != "-":
            p2 += 1

        if a != b:
            buffer.append(f"{sequence_name}\t{p1 if a != '-' else '-'}\t{a}\t{b}\t{p2 if b != '-' else '-'}\n")
            line_mismatches += 1
            if a == "-" or b == "-":
                line_gaps += 1

        if a in ("A", "C", "G", "T") and b in ("A", "C", "G", "T"):
            output += "".join(buffer)
            line_len += line_len_curr
            line_len_curr = 0
            buffer.clear()

    return output, line_len, line_mismatches, line_gaps, found_start, buffer, p1, p2


def parse(input_file, output_file, sequence_name: str, offset_seq1 = 0, offset_seq2 = 0):
    """
    Process an alignment file that is structured in blocks like:
                       10        20        30        40        50
    141486 CATAATTTTATTAGGTAGGGTCTGTGTGTCTACCACCCAATATTCTTTTG
           :::::::::::::::::::::::::::::::::::::::: :::::::::
    141434 CATAATTTTATTAGGTAGGGTCTGTGTGTCTACCACCCAAAATTCTTTTG
                   10        20        30        40        50

    Writes mismatch/gaps information into output_file in TSV format with columns:
        sequence name, position, reference base, alternate base, info.

    Notes
    -----
    - The function ignores leading bases until it finds the first valid base pair (A/C/G/T)
        in both sequences at the start and end.
     - Base mismatches are counted only for positions where both sequences
        contain a nucleotide; positions with gaps ('-') are excluded from
        the mismatch count.
    - Gaps ('-') are recorded separately in the output.

    Returns
    -------
    seq_len : int
        Total number of positions processed in the sequences.
    mismatches : int
        Total number of base mismatches (excluding gaps '-') detected between the two sequences.
    gaps : int
         Total number of gaps ('-') observed in unmatched positions.
    """

    with open(input_file) as infile, open(output_file, 'w') as outfile:
        outfile.write("sequence name\tPOS\tREF\tALT\tINFO\n")
        found_start = False
        buffer = []  # Temporarily store potential mismatches/gaps until they are confirmed by subsequent aligned bases,
                     # ensuring that mismatches occurring at the leading or trailing ends of the sequence are not reported.
        seq_len, mismatches, gaps = 0, 0, 0

        # read header
        line = infile.readline()
        while line.startswith("#"):
            line = infile.readline()

        line = infile.readline()

        while line.startswith("#"):
            line = infile.readline()
        infile.readline()
        # end read header

        curr1, curr2 = offset_seq1, offset_seq2
        while True:
            if line.startswith("#"):
                break

            seq1_line = (infile.readline().strip().split()) # skip upper positions
            infile.readline()  # skip match line
            seq2_line = (infile.readline().strip().split())
            infile.readline()  # skip positions lower line
            infile.readline()  # skip free line

            if seq1_line[0].startswith("#") or seq2_line[0].startswith("#"):
                break

            seq1_line = seq1_line[1]
            seq2_line = seq2_line[1]

            if not seq2_line:
                break

            line, line_length, mismatches_add, gaps_add, found_start, buffer, curr1, curr2 = check_sequences(
                sequence_name, [seq1_line, curr1], [seq2_line, curr2], found_start, buffer
            )
            outfile.write(line)

            seq_len += line_length
            mismatches += mismatches_add
            gaps += gaps_add

            line = infile.readline()
            if not line:
                break

        return seq_len, mismatches, gaps


def print_usage():
    print(
        "Usage:\n"
        "  python stretcher_parser.py <input_file> <output_file_tsv> <sequence_name> <offset_seq1> <offset_seq2>\n\n"
        "Example:\n"
        "  python script.py alignment.aln result.tsv chr1 0 0\n",
        file=sys.stderr
    )


if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("Error: Wrong number of arguments.\n", file=sys.stderr)
        print_usage()
        sys.exit(1)

    in_file = sys.argv[1]
    out_file = sys.argv[2]
    seq_name = sys.argv[3]

    try:
        offset1 = int(sys.argv[4])
        offset2 = int(sys.argv[5])
    except ValueError:
        print("Error: Offsets must be integers.\n", file=sys.stderr)
        print_usage()
        sys.exit(1)

    # Run parsing
    length, mismatches, gaps = parse(in_file, out_file, seq_name, int(offset1), int(offset2))

    # Print statistics
    print(f"Alignment statistics for sequence '{seq_name}':")
    print(f"  Total aligned positions (excluding ignored start and end): {length}")
    print(f"  Base mismatches (including gaps): {mismatches + gaps}")
    print(f"  Gaps detected: {gaps}")
    if length > 0:
        mismatch_rate = (mismatches + gaps) / length * 100
        gap_rate = gaps / length * 100
        print(f"  Mismatch rate (including gaps): {mismatch_rate:.8f}%")
        print(f"  Gap rate: {gap_rate:.8f}%")
