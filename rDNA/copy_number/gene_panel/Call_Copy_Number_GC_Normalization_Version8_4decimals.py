import re
import argparse
import numpy
import random
import math
import statistics
import multiprocessing as mp
import gzip
import os
import csv
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns


def get_args():
  parser = argparse.ArgumentParser(description="Find copy number of a feature in WGS data with pre-counted kmer multiplicities in the feature and a normalization set")
  parser.add_argument("-r1", help="A jellyfish counts fasta of kmers from the feature in the first WGS read (for paired and unpaired)")
  parser.add_argument("-r2", required = False, help="A jellyfish counts fasta of kmers from the feature in the second WGS read (if paired)")
  parser.add_argument("-nc1", help="A jellyfish counts fasta of kmers from the matched windows in the first WGS read (for paired and unpaired). Kmers have ideally been filtered for those unique to the matched window")
  parser.add_argument("-nc2", required = False, help="A jellyfish counts fasta of kmers from the matched windows in the second WGS read (for paired). Kmers have ideally been filtered for those unique to the matched window")
  parser.add_argument("-ngcn", help="A jellyfish counts file of kmers from the matched window, so counts are indicating the multiplicity of the kmers in the matched window.")
  parser.add_argument("-ID", help="A sample ID to use in the output")
  parser.add_argument("-gcn", help="A jellyfish counts file of kmers from the feature of interest, so counts are indicating the multiplicity of the kmers in the feature.")
  return parser.parse_args()

def takeSecond(elem):
    return int(elem[1])


args = get_args()
counts_r1 = args.r1
if args.r2 is not None:
    counts_r2 = args.r2
else:
    counts_r2 = False
norm_counts1 = args.nc1
if args.nc2 is not None:
    norm_counts2 = args.nc2
else:
    norm_counts2 = False
norm_gcn = args.ngcn
ID = args.ID
gcn= args.gcn

counts_list = []
kmer_list = []
kmer_dic = {}

with open(counts_r1, "r+") as set:
    while True:
        l1 = set.readline().strip()
        kmer = set.readline().strip()
        if l1 == "":
            break
        regex = re.search(r'\>([\S]+)',l1)
        count = int(regex.group(1))
        kmer_dic[kmer] = count
        #kmer_list.append([kmer,count])

if counts_r2 is not False:
    with open(counts_r2, "r+") as set:
        while True:
            l1 = set.readline().strip()
            kmer = set.readline().strip()
            if l1 == "":
                break
            regex = re.search(r'\>([\S]+)',l1)
            count = int(regex.group(1))
            if kmer in kmer_dic:
                kmer_dic[kmer] = kmer_dic[kmer] + count
            else:
                kmer_dic[kmer] = count
#kmer_list.sort(key=takeSecond, reverse=True)
#print(kmer_list[0:100])
#print(kmer_list[-100:])

with open(gcn, "r+") as fh:
    genomic_copy_number_dic = {}
    while True:
        l1 = fh.readline().strip()
        kmer = fh.readline().strip()
        if l1 == "":
            break
        regex = re.search(r'\>([\S]+)',l1)
        count = int(regex.group(1))
        genomic_copy_number_dic[kmer] = count


genomic_cn_normalized_counts_dic = {}
genomic_cn_normalized_counts_list = []
for kmer in kmer_dic:
    genomic_cn_normalized_counts_dic[kmer] = kmer_dic[kmer]/genomic_copy_number_dic[kmer]
    genomic_cn_normalized_counts_list.append(genomic_cn_normalized_counts_dic[kmer])

print(genomic_cn_normalized_counts_dic)


standard_dev = statistics.stdev(genomic_cn_normalized_counts_list)
mean = statistics.mean(genomic_cn_normalized_counts_list)

final_list = []
for item in genomic_cn_normalized_counts_list:
    count = item
    if count != 0:
        if abs(mean-count) < 3*standard_dev:
            final_list.append(item)
        else:
            print(item)


#find median
median = statistics.median_low(final_list)
print("The median genomic cn adjusted feature kmer counts")
print(median)
mean = statistics.mean(final_list)
print("The mean genomic cn adjusted feature kmer counts")
print(mean)

norm_counts_dic = Counter()
with open(norm_counts1, "r+") as norm1:
    while True:
        l1 = norm1.readline().strip()
        kmer = norm1.readline().strip()
        if l1 == "":
            break
        regex = re.search(r'\>([\S]+)',l1)
        count = int(regex.group(1))
        norm_counts_dic[kmer] += count

if norm_counts2 is not False:
    with open(norm_counts2, "r+") as norm2:
        while True:
            l1 = norm2.readline().strip()
            kmer = norm2.readline().strip()
            if l1 == "":
                break
            regex = re.search(r'\>([\S]+)',l1)
            count = int(regex.group(1))
            norm_counts_dic[kmer] += count

norm_gcn_dic = Counter()
with open(norm_gcn) as ngcn:
    while True:
        l1 = ngcn.readline().strip()
        kmer = ngcn.readline().strip()
        if l1 == "":
            break
        regex = re.search(r'\>([\S]+)',l1)
        count = int(regex.group(1))
        norm_gcn_dic[kmer] += count

norm_counts_adj = Counter()
for item in norm_counts_dic:
    print(item)
    print(norm_counts_dic[item])
    print(norm_gcn_dic[item])
    norm_counts_adj[item] = float(norm_counts_dic[item]/norm_gcn_dic[item])


print("The mean normalization set count")
print(statistics.mean(list(norm_counts_adj.values())))
norm_mean = statistics.mean(list(norm_counts_adj.values()))
diploid_norm_mean = norm_mean/2

norm_median = statistics.median(list(norm_counts_adj.values()))
print("The median normalization set count")
print(norm_median)
diploid_norm_median = norm_median/2

ribosome_copy_number = median/norm_mean
print("Median Copy Number" + "\n")
print(ribosome_copy_number)

mean_copy_number = mean/norm_mean
print("Mean Copy Number:" + "\n")
print(mean_copy_number)

diploid_ribosome_copy_number = median/diploid_norm_mean
print("Diploid Median Copy Number" + "\n")
print(diploid_ribosome_copy_number)

diploid_mean_copy_number = mean/diploid_norm_mean
print("Diploid Mean Copy Number:" + "\n")
print(diploid_mean_copy_number)

exists = os.path.isfile("./Copy_Numbers.tsv")

summary = open("Copy_Numbers.tsv","a+")

if exists == False:
    summary.write("Sample" + "\t" + "Median_Haploid_Copy_Number" + "\t" + "Median_Diploid_Copy_Numer" + "\t" + "Mean_Haploid_Copy_Number" + "\t" + "Mean_Diploid_Copy_Number" + "\n")
    summary.write(ID + "\t" +str(round(ribosome_copy_number,4)) + "\t" + str(round(diploid_ribosome_copy_number,4)) + "\t" + str(round(mean_copy_number,4)) + "\t" + str(round(diploid_mean_copy_number,4))+"\n")
else:
    summary.write(ID + "\t" +str(round(ribosome_copy_number,4)) + "\t" + str(round(diploid_ribosome_copy_number,4)) + "\t" + str(round(mean_copy_number,4)) + "\t" + str(round(diploid_mean_copy_number,4))+"\n")
summary.close()



#plt.hist(norm_counts_dic.values())
sns.displot(norm_counts_dic.values())
plt.savefig("adjusted_matched_window_counts.png")
plt.show()

raw_feature_counts_list = [int(i) for i in kmer_dic.values()]
outfil = open("raw_feature_counts.txt", "w+")
for item in raw_feature_counts_list:
    outfil.write(str(item) + "\n")
outfil.close()
sns.displot(raw_feature_counts_list)
#plt.hist(kmer_dic.values(), range=(1,max(raw_feature_counts_list)+10))
plt.savefig("raw_feature_counts.png")
plt.show()

genomic_cn_normalized_counts_list = [int(i) for i in genomic_cn_normalized_counts_list]
outfil = open("adjusted_feature_counts.txt", "w+")
for item in genomic_cn_normalized_counts_list:
    outfil.write(str(item) + "\n")
outfil.close()
sns.displot(genomic_cn_normalized_counts_list)
#plt.hist(genomic_cn_normalized_counts_list)
plt.savefig("adjusted_feature_counts.png")
plt.show()
