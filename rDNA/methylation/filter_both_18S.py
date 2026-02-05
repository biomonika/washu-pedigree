import sys
import re

def overlap(curr_region, wanted_region):
	if curr_region[1] < wanted_region[0] or curr_region[0] > wanted_region[1]:
		return 0
	overlap = min(curr_region[1], wanted_region[1]) - max(curr_region[0] - wanted_region[0])
	return overlap / (wanted_region[1] - wanted_region[0])

def subregion_ani(t_region, wanted_region, cigar):
	cigar_pattern = r'(\d+)([MIDNSHP=X])'
	cigar_tuples = re.findall(cigar_pattern, cigar)
	ops, lengths = zip(*cigar_tuples)
	curr_pos = t_region[0]
	matching_bases = 0
	for i in range(len(ops)):
		op, length = ops[i], lengths[i]
		if curr_pos >= wanted_region[0] and curr_pos <= wanted_region[1]:
			if op in ['M', '==', 'X']:
				matching_bases += length
		if op in ['M', '=', 'X', 'D', 'N']:
			curr_pos += length
	return matching_bases / (wanted_region[1] - wanted_region[0])

sam_no_header_f = open(sys.argv[1])
paf_f = open(sys.argv[2])
output_type = sys.argv[3]
if len(sys.argv) == 5:
	pi_cutoff = int(sys.argv[4])/100
	aln_cutoff = int(sys.argv[4])/100
elif len(sys.argv) == 6:
	pi_cutoff = int(sys.argv[4])/100
	aln_cutoff = int(sys.argv[5])/100
else:
	pi_cutoff = 0.9
	aln_cutoff = 0.9

want_subregion=False
verbose=False
if len(sys.argv) > 6:
	if '.bed' in sys.argv[5]:
		want_subregion=True
		bed_f = open(sys.argv[5])
		for line in bed_f:
			line_split = line.strip().split('\t')
			subregion = (int(line_split[1]), int(line_split[2]))
		if len(sys.argv) > 7:
			verbose=True
	else:
		verbose=True

regex_str = 'de:f:((0|\d+)\.?\d*)'
for paf_line, sam_line in zip(paf_f, sam_no_header_f):
	if verbose: print(sam_line.strip())
	paf_line_split = paf_line.strip().split('\t')
	sam_line_split = sam_line.strip().split('\t')
	t_start, t_end = int(paf_line_split[7]), int(paf_line_split[8])
	cigar = re.findall('cg:Z:(\w+)', paf_line)[0]
	if want_subregion:
		if overlap((t_start, t_end), subregion) >= pi_cutoff:
			if subregion_ani((t_start, t_end), subregion, cigar) >= pi_cutoff:
				if output_type == 'sam':
					print(sam_line.strip())
				elif output_type == 'paf':
					print(paf_line.strip())
	else:
		q_len, matching_bases, aln_block = int(paf_line_split[6]), int(paf_line_split[9]), int(paf_line_split[10])
		de_flag = re.findall(regex_str, sam_line)[0][0]
		ani = 1 - float(de_flag)
		if verbose: print(paf_line.strip())
		if verbose: print(aln_block, aln_cutoff*q_len, ani)
		if aln_block >= aln_cutoff * q_len and ani >= pi_cutoff:
			if verbose: print('keep')
			if output_type == 'sam':
				print(sam_line.strip())
			elif output_type == 'paf':
				print(paf_line.strip())
		if verbose: print('-----')
