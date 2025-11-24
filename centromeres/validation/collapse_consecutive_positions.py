import sys

def flush_buffer(buf):
    """
    Consolidate buffered elements into a single output line.

    Parameters
    ----------
    buf : list of lists
        Each element is an array: [seq, pos_ref, ref, alt, pos_alt] representing
        a mismatch or gap.

    Returns
    -------
    str
        A single string combining the buffer information into one TSV-formatted line.
    """

    if not buf:
        return ""

    seq = buf[0][0]
    pos = buf[0][1]
    ref = "".join(row[2] for row in buf if row[2] != "-")
    alt = "".join(row[3] for row in buf if row[3] != "-")
    info = buf[0][4]

    return f"{seq}\t{pos}\t{ref}\t{alt}\t{info}\n"


def is_consecutive(prev_ref, prev_alt, ref, alt):
    """
    Determine whether the current position is consecutive with the previous position from the input TSV file.

    Consecutive positions are defined as positions that increment by one for either
    or both sequences, while properly handling gaps ('-'):
        - Both reference and alternate positions advance by 1.
        - Reference-only advance: alternate remains a gap.
        - Alternate-only advance: reference remains a gap.

    Parameters
    ----------
    prev_ref, prev_alt : str
        Previous reference and alternate positions (may be '-' for gaps)
    ref, alt : str
        Current reference and alternate positions (may be '-' for gaps)

    Returns
    -------
    bool
        True if the current pair is consecutive, False otherwise.
    """

    # both advance
    if prev_ref != "-" and ref != "-" and prev_alt != "-" and alt != "-":
        return int(ref) == int(prev_ref) + 1 and int(alt) == int(prev_alt) + 1
    # reference-only advance (alt must stay '-')
    if prev_ref != "-" and ref != "-" and alt == "-" and prev_alt == "-":
        return int(ref) == int(prev_ref) + 1
    # alt-only advance (ref must stay '-')
    if prev_alt != "-" and alt != "-" and ref == "-" and prev_ref == "-":
        return int(alt) == int(prev_alt) + 1
    return False


def reformat(input_file, output_file):
    """
    Reformat a TSV file (output from stretcher_parser.py) with mismatches and gaps, collapsing consecutive positions
    into a single line.

    Parameters
    ----------
    input_file : str
        Path to the input TSV file. Expected columns: seq, pos_ref, ref, alt, pos_alt
    output_file : str
        Path to write the reformatted TSV file.

    Notes
    -----
    - The output file will contain columns with the information in following order:
        - sequence name
        - reference position
        - reference nucleotide or gap ('-')
        - alternate position
        - alternate nucleotide or gap ('-')
    """

    with open(input_file, 'r') as in_fh, open(output_file, 'w') as out_fh:
        current_seq = None
        buffer = []

        # Write header
        out_fh.write("seq\tpos\tref\talt\tinfo\n")

        for line in in_fh:
            row = line.strip().split("\t")
            if not row or len(row) < 5:
                continue

            seq, pos_ref, ref, alt, pos_alt = row

            # Initialize sequence tracking
            if not current_seq:
                current_seq = seq
                continue

            # Start buffering if empty
            if len(buffer) == 0:
                buffer.append(row)
                continue

            pos_ref_prev, pos_alt_prev = buffer[-1][1], buffer[-1][4]

            if is_consecutive(pos_ref_prev, pos_alt_prev, pos_ref, pos_alt):
                buffer.append(row)
            else:
                # Write out buffered consecutive block
                out_fh.write(flush_buffer(buffer))
                buffer = [row]
                current_seq = seq

        # Write remaining buffer
        out_fh.write(flush_buffer(buffer))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python reformat_tsv.py <input_file_tsv> <output_file_tsv>", file=sys.stderr)
        sys.exit(1)

    in_file = sys.argv[1]
    out_file = sys.argv[2]
    reformat(in_file, out_file)
