#!/usr/bin/env python3
import sys

"""
- reads in bed file containing regions of one or two transmissions
- sorts the regions
- merges start and end of the regions if a certain threshold is not met (1 Mb)
"""

def main(bedfile: str, outfile: str) -> None:

    chrom = ""
    haplo = ""

    intervals = []
    with open(bedfile) as bedf:
        for line in bedf:
            parts = line.strip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            intervals.append((start, end))

    intervals.sort()
    newintervals = []

    newstart = intervals[0][0]
    newend = intervals[0][1]

    for i in range(1, len(intervals)):
        start, end = intervals[i]

        if newend > end:
            newend = end
        else:
            newintervals.append((newstart, newend))
            newstart = start
            newend = end

    newintervals.append((newstart, newend))

    with open(outfile, "w") as outf:
        for s, e in newintervals:
            outf.write(chrom + "\t" + haplo + "\t" + str(s) + "\t" + str(e) + "\n")


if __name__ == "__main__":
    bedfile = sys.argv[1]
    if len(sys.argv) > 2:
        outfile = sys.argv[2]
    else:
        outfile = ".".join(bedfile.split(".")[:-1]) + "_merged.bed"

    main(bedfile, outfile)
