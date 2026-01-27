#!/usr/bin/env python3
import sys
from cyvcf2 import VCF


def main(phasedfile: str, outfile: str) -> None:
	with open(outfile, "w") as outf:
		outf.write(
				"chr\t#switches\tlist of switches\n"
		)

		chromnames = []
		chromosomes = VCF(phasedfile).seqnames

		for chrom in chromosomes:
			print(chrom)
			ps_to_calls = {}
			for variant in VCF(phasedfile):
				if variant.CHROM != chrom:
					continue
	
				ps = variant.format("PS")[0][0]
				if ps not in ps_to_calls.keys():
					ps_to_calls[ps] = []
				ps_to_calls[ps].append(variant)
	
			for ps, varlist in ps_to_calls.items():
				cur = varlist[0].genotypes[0][0]
				curpos = varlist[0].POS
				switches_temp = [curpos]
	
				for var in varlist[1:]:
					if var.genotypes[0][0] != cur:  # phase switch
						if len(switches_temp) < 2:
							switches_temp.append(var.POS)
						else: # if the distance between switches is less than 100 kb
							if var.POS - switches_temp[-1] < 100000:
								switches_temp.pop()
							else:
								switches_temp.append(var.POS)
				cur = var.genotypes[0][0]
	
				switches_temp.append(varlist[-1].POS)
				print(switches_temp, ps)
				if len(switches_temp) > 2:
					chromnames.append(chrom)
	
				outf.write(chrom+"\t"+str(len(switches_temp))+"\t"+",".join(str(p) for p in switches_temp)+"\n")
		print(chromnames)	

if __name__ == "__main__":
	phasedfile = sys.argv[1]
	outfile = sys.argv[2]
	main(phasedfile, outfile)
