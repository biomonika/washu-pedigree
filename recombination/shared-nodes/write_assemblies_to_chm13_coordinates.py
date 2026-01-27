#!/usr/bin/env python3
import sys

"""
- takes a bed file with the regions in the assembly where three-generational nodes map to
- takes a fragment ID from the assembly (e.g. haplotype1-0000003) as input and a chromosome (chr01) determined by the ragtag agp file
- writes an output .bed file containing the mapped regions
"""

def chrom_list_1_to_22():
    return [f"chr0{i}" for i in range(1, 10)] + [f"chr{i}" for i in range(10, 23)]


def normalize_chrom(chrom: str) -> str:
    
    if len(chrom) == 4:  # "chr1".."chr9"
        return chrom.replace("chr", "chr0")
    return chrom


def parse_ragtag_agp(agpfile: str):
    """
    Returns: dict chrom -> list of (fragment_id, offset)
      - only lines starting with "chr"
      - skip rows where parts[6] == "scaffold"
      - chrom name is parts[0].split('_')[0]
      - offset is int(parts[1]); if offset == 1, use 0
    """
    chrom_to_fragments = {}

    with open(agpfile) as agp:
        for line in agp:
            if not line.startswith("chr"):
                continue

            parts = line.strip().split("\t")
            chrom = normalize_chrom(parts[0].split("_")[0]) #remove trailing _RagTag suffix from the chrom id

            if chrom not in chrom_to_fragments:
                chrom_to_fragments[chrom] = []

            if parts[6] == "scaffold":
                continue

            frag = parts[5]
            offset = int(parts[1])
            if offset == 1:
                offset = 0
            chrom_to_fragments[chrom].append((frag, offset))

    return chrom_to_fragments


def main(bedfile: str, agpfile: str, outpath: str) -> None:
    chroms = chrom_list_1_to_22()
    chrom_to_fragments = parse_ragtag_agp(agpfile)

    for chrid in chroms:
        outfile_scaf = outpath + "mapped_boundaries_" + chrid + "_scaffold.bed"

        with open(outfile_scaf, "w") as outf_scaf:
            if chrid not in chrom_to_fragments:
                outf_scaf.write("")
                continue

            for fragmentid, offset in chrom_to_fragments[chrid]:
                with open(bedfile) as bed:
                    for line in bed:
                        parts = line.strip().split("\t")
                        fr = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])

                        if fr == fragmentid:
                            newstart = offset + start
                            newend = offset + end
                            outf_scaf.write(
                                chrid
                                + "\t"
                                + fr
                                + "\t"
                                + str(newstart)
                                + "\t"
                                + str(newend)
                                + "\n"
                            )


if __name__ == "__main__":
    bedfile = sys.argv[1]
    agpfile = sys.argv[2]
    outpath = sys.argv[3]
    main(bedfile, agpfile, outpath)
