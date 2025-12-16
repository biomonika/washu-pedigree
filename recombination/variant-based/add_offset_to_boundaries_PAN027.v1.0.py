#!/usr/bin/env python3
import sys


def normalize_chr(name: str) -> str:
    # "chr1" -> "chr01", "chr10" -> "chr10", "chrX" stays "chrX"
    if name != "chrX" and int(name.split("chr")[1]) < 10:
        return "chr0" + name.split("chr")[1]
    return name


def chromosomes_1_to_22() -> list[str]:
    return [f"chr0{i}" for i in range(1, 10)] + [f"chr{i}" for i in range(10, 23)]


def main(agpfile: str, filepath: str) -> None:
    allchromosomes = chromosomes_1_to_22()
    print(allchromosomes)

    for chrom in allchromosomes:
        infile = f"{filepath}/recomb_boundaries_{chrom}.bed"
        phaseset_infile = f"{filepath}/phaseset_boundaries_{chrom}.bed"
        scaffolded_boundaryfile = f"{filepath}/recomb_boundaries_{chrom}_scaffold.bed"
        scaffolded_phasesetfile = f"{filepath}/phaseset_boundaries_{chrom}_scaffold.bed"
        contigfile = f"{filepath}/contig_boundaries_{chrom}.bed"

        boundarylines = []
        phasesetlines = []
        contiglines = {}

        contigstart, contigend = -1, -1

        # read offset from ragtag_scaffold.agp file
        offset = 0
        linecontig = None

        with open(agpfile) as agp:
            for line in agp:
                parts = line.strip().split("\t")
                if parts[0] != chrom:
                    continue

                linecontig = normalize_chr(parts[5].split(".")[1])
                if linecontig != chrom:
                    continue

                lineoffset = int(parts[1])
                offset = lineoffset - 1  # .agp file intervals seem to be 1-indexed
                contigstart = offset
                contigend = int(parts[2])
                contiglines[chrom] = (
                    parts[0]
                    + "\t"
                    + linecontig
                    + "\t"
                    + str(contigstart)
                    + "\t"
                    + str(contigend)
                    + "\n"
                )

        # read boundaries from infile, add offset, and write scaffolded output
        with open(infile) as inf:
            for line in inf:
                parts = line.strip().split("\t")
                start = int(parts[2])
                end = int(parts[3])
                newstart = start + offset
                newend = end + offset
                boundarylines.append(
                    parts[0]
                    + "\t"
                    + linecontig
                    + "\t"
                    + str(newstart)
                    + "\t"
                    + str(newend)
                    + "\t"
                    + parts[-1]
                    + "\n"
                )

        with open(phaseset_infile) as ps_inf:
            for line in ps_inf:
                parts = line.strip().split("\t")
                start = int(parts[1])
                end = int(parts[2])
                newstart = start + offset
                newend = end + offset
                phasesetlines.append(
                    parts[0]
                    + "\t"
                    + str(newstart)
                    + "\t"
                    + str(newend)
                    + "\t"
                    + parts[-1]
                    + "\n"
                )

        with open(scaffolded_boundaryfile, "w") as scaf:
            for line in boundarylines:
                scaf.write(line)

        with open(scaffolded_phasesetfile, "w") as ps:
            for line in phasesetlines:
                ps.write(line)

        with open(contigfile, "w") as contigf:
            contigf.write(contiglines[chrom])


if __name__ == "__main__":
    agpfile = sys.argv[1]
    filepath = sys.argv[2]
    main(agpfile, filepath)
