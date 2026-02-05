import sys

in_f=open(sys.argv[1])
grep_f=open(sys.argv[2])

chimeras=set()
for line in grep_f:
	chimeras.add(line.strip())

for line in in_f:
	if line.split('\t')[0] not in chimeras:
		print(line.strip())
