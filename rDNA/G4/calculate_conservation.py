#!/usr/bin/env python3
# calculate conservation of an MSA
# how many basepairs are matching the consensus represented by the first sequence?

from Bio import AlignIO
import sys

if len(sys.argv) != 4:
    sys.exit(f"Usage: {sys.argv[0]} <msa.fa> <output.bed> <seq_name>")

msa_file = sys.argv[1]
output_file = sys.argv[2]
seq_name = sys.argv[3]

aln = AlignIO.read(msa_file, "fasta")
ref = aln[0].seq

with open(output_file, "w") as out:
    for pos in range(aln.get_alignment_length()):
        ref_base = ref[pos].upper()

        # Skip unknown bases or gaps in the consensus
        if ref_base in ("N", "-"):
            continue

        score = 0
        for rec in aln[1:]:  # skipping the consensus as the first sequence
            base = rec.seq[pos].upper()

            # Ignore gaps and unknown bases
            if base in ("-", "N"):
                continue

            if base == ref_base:
                score += 1

        # BED-like output: seq_name, 0-based position, end, ref base, conservation score
        out.write(
            f"{seq_name}\t{pos}\t{pos+1}\t{ref_base}\t{score}\n"
        )