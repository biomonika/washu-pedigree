#!/usr/bin/env python3
import sys


def normalize_chrom(component: str) -> str:
    # component like "chr1", "chr10", "chrX"
    if component != "chrX" and int(component.split("chr")[1]) < 10:
        return "chr0" + component.split("chr")[1]
    return component


def main(infile: str, outpath_prefix: str, phasesetpath_prefix: str) -> None:
    # Input file must not contain spaces after commas in the list of positions
    chroms = []
    chrom_to_lines = {}
    chrom_to_phasesets = {}

    chrom = None
    component = None

    with open(infile) as inf:
        for line in inf:
            if line.startswith("P"):  # e.g. 'PAN027.chr1.maternal'
                component = line.strip().split(".")[1]  # e.g. "chr1"
                chrom = normalize_chrom(component)

                chroms.append(chrom)
                if chrom not in chrom_to_lines:
						chrom_to_lines[chrom] = []
					if chrom not in chrom_to_phasesets:
						chrom_to_phasesets[chrom] = []
               continue

            # positions line, e.g. "[138030,8113343,...]"
            vals = (
                line.strip()
                .split(" ")[0]
                .strip("[")
                .strip("]")
                .split(",")
            )

            chrom_to_phasesets[chrom].append(vals)

            for i in range(len(vals) - 1):
                col = "0" if i % 2 == 0 else "1"

                # If the value pair is the first, add the interval from 0 to first boundary before
                if len(chrom_to_lines[chrom]) == 0:
                    chrom_to_lines[chrom].append(
                        chrom + "\t" + component + "\t" + "0" + "\t" + str(vals[i]) + "\t" + "2" + "\n"
                    )

                chrom_to_lines[chrom].append(
                    chrom + "\t" + component + "\t" + str(vals[i]) + "\t" + str(vals[i + 1]) + "\t" + col + "\n"
                )

    print(chroms)

    for c in sorted(chroms):
        outfile = outpath_prefix + "_" + c + ".bed"
        with open(outfile, "w") as outf:
            for line in chrom_to_lines[c]:
                outf.write(line)

        phaseset_outfile = phasesetpath_prefix + "_" + c + ".bed"
        with open(phaseset_outfile, "w") as phaseout:
            vals = chrom_to_phasesets[c]
            if len(vals) > 0:
                for v in vals:
                    phaseout.write(c + "\t" + str(v[-1]) + "\t" + str(v[-1]) + "\t" + "v" + "\n")


if __name__ == "__main__":
    infile = sys.argv[1]
    outpath = sys.argv[2]
    phasesetpath = sys.argv[3]
    main(infile, outpath, phasesetpath)
