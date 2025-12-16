import sys

"""
- takes the file with three-generational nodes from the paternal (maternal) side and the paf file from its mapping to the assemblies
- writes the regions where the contigs map into bed file
- outputs bed file with the mapped regions in the assembly coordinates 

- currently only takes contigs with one unique match 
"""

contigfile = sys.argv[1] 
paffile = sys.argv[2] 
outfile = sys.argv[3]

contigs = []
with open(contigfile) as inf:
	for i,l in enumerate(inf):
		contigs.append(l.strip())

lowqual = []
mapped_contigs = {}
with open(paffile) as paf:
	for i,l in enumerate(paf):
		parts = l.strip().split('\t')
		contig = parts[0]
		target = parts[5]
		start = parts[7]
		end = parts[8]
		mapq = int(parts[11])
		if mapq != 60:
			lowqual.append(lowqual)
		if mapq > 0:
			if contig not in mapped_contigs:
				mapped_contigs[contig] = [] 
			mapped_contigs[contig].append((target, start, end))
		
multiple= []		
with open(outfile, 'w') as outf:		
	for contig,intervals in mapped_contigs.items():
		if len(intervals)==1:
			target = intervals[0][0]
			start = intervals[0][1]
			end = intervals[0][2]
			outf.write(target+'\t'+start+'\t'+end+'\t'+contig)
			outf.write('\n')
		else:
			multiple.append(contig)
			


