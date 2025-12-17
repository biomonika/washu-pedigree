import re
import argparse
import numpy
import random
import math
import statistics
import multiprocessing as mp

def get_args():
  parser = argparse.ArgumentParser(description="Filter the single copy set for the samples being assessed. Some kmers may not exist/ no longer be single copy")
  parser.add_argument("-fa")
  parser.add_argument("-ID")
  parser.add_argument("-read")
  parser.add_argument("-k")
  return parser.parse_args()

args = get_args()
fa = args.fa
ID = args.ID
read = args.read
k = args.k

def takeSecond(elem):
    return int(elem[1])



with open(fa, "r+") as set:
    kmer_list = []
    while True:
        l1 = set.readline().strip()
        kmer = set.readline().strip()
        if l1 == "":
            break
        regex = re.search(r'\>([\S]+)',l1)
        count = regex.group(1)
        kmer_list.append([kmer,count])

print(kmer_list[0:10])
kmer_list.sort(key=takeSecond, reverse=True)
#sorted_list_of_counts = aggregated_counts.sort(key=takeThird)
#sorted_list_of_counts = list_of_counts.sort(reverse=True,key=lambda x: x[2])

sum = 0
kmers_counted = 0
stdev_list = []
for item in kmer_list:
    if int(item[1]) != 0:
        sum += int(item[1])
        kmers_counted += 1
        stdev_list.append(int(item[1]))
mean_scs = sum/kmers_counted
standard_dev_scs = statistics.stdev(stdev_list)

print(mean_scs)
print(standard_dev_scs)

final_list = []
for item in kmer_list:
    count = int(item[1])
    if count != 0:
        if abs(mean_scs-count) < 3*standard_dev_scs:
            final_list.append(item)
        else:
            print(item)

print(mean_scs)
print(standard_dev_scs)
comparison_list = []
for item in final_list:
    count = int(item[1])
    comparison_list.append(count)
print(statistics.mean(comparison_list))
print(statistics.stdev(comparison_list))


string = "matched_gc_" + ID + "_" + "k" + k + "_" + "filtered_" + read + ".fa"
new_fasta = open(string, "w+")
for number in range(len(final_list)):
    new_fasta.write(">")
    new_fasta.write(str(final_list[number][1]))
    new_fasta.write("\n")
    new_fasta.write(str(final_list[number][0]))
    new_fasta.write("\n")
