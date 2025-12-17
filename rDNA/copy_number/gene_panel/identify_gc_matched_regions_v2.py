import re
import argparse
import numpy
import random
import math
from collections import Counter

def get_args():
  parser = argparse.ArgumentParser(description="Find Regions of the Genome wihin 1% GC content of target sequence. Also avoids rDNA (project specific criteria)")
  parser.add_argument("-rbed", required=False)
  parser.add_argument("-assembly")
  parser.add_argument("-roi")
  parser.add_argument("-chr", help="A list of chromosome IDs for which to extract gc matched sequence", required= False)
  parser.add_argument("-w_size", help="The size of the normalization bins to tile the genome with", default = 2000)
  return parser.parse_args()

args = get_args()
if args.rbed is not None:
    bed = args.rbed
else:
    bed = False
assembly = args.assembly
roi = args.roi
if args.chr is not None:
    chrids = args.chr
else:
    chrids = False
win_size = int(args.w_size)



def gc_content(sequence):
    AT = 0
    GC = 0
    No_Ns = True
    for character in sequence:
        if character == "A" or character == "a":
            AT += 1
        elif character == "T" or character == "t":
            AT += 1
        elif character == "G" or character == "g":
            GC += 1
        elif character == "C" or character == "c":
            GC += 1
        else:
            print("Non ATGC character detected")
            No_Ns = False
            break
    if No_Ns == True:
        total = GC+AT
        ratio = GC/total
        return(round(ratio,3))
    else:
        return("Null")

with open(roi,"r+") as ref:
    while True:
        seqname = ref.readline().strip()
        if seqname == "":
            break
        seq = ref.readline().strip()
        #print(seq)
        target_gc = gc_content(seq)


#Avoid pulling from the rDNA regions of the acrocentrics
#Put these regions in a dictionary, if on acrocentric check if dic[chr] is not in range(rDNA positions)
rdna_region_list = []
if bed != False:
    with open(bed) as rdna:
        for line in rdna:
            content = line.strip()
            bed_fields = re.search(r'^([\S]+)\t([\S]+)\t([\S]+)',content)
            chrm = bed_fields.group(1)
            bed_entry = (bed_fields.group(1),int(bed_fields.group(2)),int(bed_fields.group(3)))
            rdna_region_list.append(bed_entry)

chrID_list = []
if chrids != False:
    with open(chrids) as fil:
        for line in fil:
            chrID_list.append(line.strip())

valid_window_list = []
with open(assembly, "r+") as ref:
    while True:
        seqname = ref.readline().strip()
        if seqname == "":
            break
        seq = ref.readline().strip()
        seqID = re.search(r'\>([\S]+)',seqname)
        if seqID.group(1) is None:
            print(seqname)
        chrmID = seqID.group(1)
        if chrids != False:
            if chrmID not in chrID_list:
                if ">" + chrmID not in chrID_list:
                    continue
        for num in range(int(math.floor(len(seq)/win_size))):
            window = seq[(0+win_size*num):(win_size-1+win_size*num)]
            GC_ratio = gc_content(window)
            if GC_ratio != "Null":
                if abs(target_gc - GC_ratio) <= 0.01:
                    start = win_size*num
                    end = win_size-1+win_size*num
                    in_blacklist = False
                    if bed != False:
                        for n in range(len(rdna_region_list)):
                            #check if the start and end of the window falls within any rdna regions
                            #if so, move on to a new window using break
                            if rdna_region_list[n][0] == chrmID:
                                if start > rdna_region_list[n][1] and start < rdna_region_list[n][2]:
                                    print("this window is in the blacklisted region!")
                                    print(chrmID)
                                    print(start)
                                    print(end)
                                    print(rdna_region_list[n][1])
                                    print(rdna_region_list[n][2])
                                    in_blacklist = True
                                    break
                                if end > rdna_region_list[n][1] and end < rdna_region_list[n][2]:
                                    in_blacklist = True
                                    break
                    if in_blacklist == False:
                        valid_window_list.append([chrmID,0+win_size*num,win_size-1+win_size*num])

outwin = open("matched_windows.bed","w+")
#outwin.write("chr"+"start"+"\t"+"end")
for n in range(len(valid_window_list)):
    outwin.write(valid_window_list[n][0]+"\t"+str(valid_window_list[n][1])+"\t"+str(valid_window_list[n][2])+"\n")
outwin.close()

#select a random subset of regions from this list to use for normalization
window_subset = []

random.seed(4)
outfile = open("matched_windows_subset.bed","w+")
if chrids != False:
    for chrm in chrID_list:
        ranges_added = 0
        current_chr = []
        for v in valid_window_list:
            if v[0] == chrm:
                current_chr.append(v)
        if len(current_chr) > 30:
            while ranges_added < 30:
                selection = random.choice(current_chr)
                window_subset.append(selection)
                ranges_added += 1
        else:
            for entry in current_chr:
                window_subset.append(entry)
                ranges_added += 1
else:
    ranges_added = 0
    max_ranges = 50
    first=True
    if len(valid_window_list) < 50:
        max_ranges = len(valid_window_list)
    while ranges_added < 50:
        selection = random.choice(valid_window_list)
        if first == True:
            prev_chr = selection[0]
            first = False
        if selection[0] != prev_chr:
            window_subset.append(selection)
            ranges_added += 1


for n in range(len(window_subset)):
    outfile.write(window_subset[n][0]+"\t"+str(window_subset[n][1])+"\t"+str(window_subset[n][2])+"\n")
outfile.close()
