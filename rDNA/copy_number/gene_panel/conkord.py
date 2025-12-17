import re
import os
import argparse
import multiprocessing
import subprocess

def get_args():
  parser = argparse.ArgumentParser(description="Find the copy number of a feature in the genome")
  parser.add_argument("--no_uniq", action='store_true', default = False)
  #parser.add_argument("-n", "--name", help="The name of the feature(s)", default = "placeholder")
  parser.add_argument("-k", nargs='+', help="Kmer size(s). Default of 31", default = 31)
  parser.add_argument("-f", nargs='+', help="A consensus sequence or single occurence for your feature(s) of interest", required = True)
  parser.add_argument("-bed", nargs='+', help="The bed file(s) with coordinates for your feature in the reference genome. Required if no_uniq is not set.")
  parser.add_argument("-r", help="The directory containing reads to be CN assessed", required = True)
  parser.add_argument("-g", nargs='+', help="The reference genome to use", required = True)
  parser.add_argument("--gzip", action='store_true', default=False)
  parser.add_argument("-t", default = 10, help = "The number of threads to use for Jellyfish kmer counting")
  parser.add_argument("-w_size", default=2000, help = "The size of genomic bins for G/C normalization")
  parser.add_argument("--cluster", action='store_true', default=False)
  return parser.parse_args()

args = get_args()
#name = args.n
k = args.k
feature_ref = sorted(args.f)
reads = args.r
genome = sorted(args.g)
if args.no_uniq == True:
    no_uniq = True
else:
    no_uniq = False
if args.bed is not None:
    fbed = sorted(args.bed)
else:
    fbed = ["NA"]
    if no_uniq == False:
        print("You must provide a bedfile with feature coordinates unless the no_uniq paramter is set! This enables finding unique feature kmers.")
        fail
#fbed = args.bed
if args.gzip ==True:
    gzip = True
else:
    gzip = False
threads = args.t
w_size = args.w_size
cluster = args.cluster

if args.bed is not None:
    for bed in fbed:
        with open(bed, "r+") as coords:
            for line in coords:
                content = line.strip()
                bed_fields = re.search(r'^([\S]+)\t([\S]+)\t([\S]+)',content)
                start = int(bed_fields.group(2))
                end = int(bed_fields.group(3))
                span = end - start
                for ksize in k:
                    if span < int(ksize):
                        print("Warning: K value is larger than at least one bed entry for your feature")


os.system('mkdir jellyfish_files results matched_windows feature')


#Get the sample IDs from the read names, propogate into config for rest of pipeline
command = 'ls ' + str(reads) + '| grep ".fq\|.fastq\|.fasta\|.fa"'
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,universal_newlines=True)
#print(p.communicate())
files = []
while True:
  line = p.stdout.readline()
  if line != "":
      files.append(line)
  if not line:
    break

#pout = p.communicate()[0]
#print(pout.split("''"))

#all_files = os.popen('ls' + string(reads) + '| grep .fq\|.fastq\|.fasta\|.fa').readlines()
#all_files = all_files.strip()


#Identify ID names from file names
id_list = []
for item in files:
    line = item.strip()
    print(line)
    regex = re.search(r'([\S,_]+)_[0-9]+\.(fa|fasta|fq|fastq|fastq.gz|fasta.gz|fa.gz|fq.gz)',line)
    id = regex.group(1)
    id_list.append(id)
#take only unique file names
id_list = set(id_list)
id_list = sorted(list(id_list))


#Do the same but for feature IDs
print(feature_ref)
print(args.f)
#if no_uniq == False:
feature_list = []
for item in feature_ref:
    line = item.strip()
    print(line)
    regex = re.search(r'/*([a-z,A-Z,0-9,-,.,_]+)\.(fa|fasta)',line)
    feature = regex.group(1)
    feature_list.append(feature)
#take only unique file names
feature_list = set(feature_list)
feature_list = sorted(list(feature_list))

print(feature_list)


# #Do the same but for feature IDs
# if no_uniq == False:
#     feature_list = []
#     for item in fbed:
#         line = item.strip()
#         print(line)
#         regex = re.search(r'/*([a-z,A-Z,0-9,-,.,_]+)\.(bed)',line)
#         feature = regex.group(1)
#         feature_list.append(feature)
#     #take only unique file names
#     feature_list = set(feature_list)
#     feature_list = list(feature_list)


#create a config file with contents adjusted by the argparse options
with open("config.yml","a+") as fil:
    fil.write("#config file for CONKORD snakemake pipeline using user parameters" + "\n")
    fil.write("FPATH:" + "\n")
    for feature in feature_ref:
        fil.write("  - \"" + feature + "\"" + "\n")
    fil.write("K:" + "\n")
    for item in k:
        fil.write("  - \"" + str(item) + "\"" + "\n")
    fil.write("ID:" + "\n")
    for id in id_list:
        fil.write("  - \"" + id + "\"" + "\n")
    fil.write("GENOME:" + "\n")
    for path in genome:
        fil.write("  - \"" + path + "\"" + "\n")
    fil.write("UNIQUE:" + "\n")
    if no_uniq == False:
        fil.write("  - \"" + "unique" + "\"" + "\n")
    if no_uniq == True:
        fil.write("  - \"" + "no_unique" + "\"" + "\n")
    fil.write("BED:" + "\n")
    for bed in fbed:
        fil.write("  - \"" + bed + "\"" + "\n")
    #if no_uniq == False:
    fil.write("FEATURE:" + "\n")
    for feature in feature_list:
        fil.write("  - \"" + feature + "\"" + "\n")
    fil.write("SEQ_DIR:" + "\n")
    fil.write("  - \"" + reads + "\"" + "\n")
    fil.write("GZIP:" + "\n")
    if gzip == True:
        fil.write("  - \"" + "True" + "\"" + "\n")
    if gzip == False:
        fil.write("  - \"" + "False" + "\"" +"\n")
    fil.write("THREADS:"+ "\n")
    fil.write("  - \"" + str(threads) + "\"" +"\n")
    fil.write("WINDOW_SIZE:"+ "\n")
    fil.write("  - \"" + str(w_size) + "\"" +"\n")
    fil.write("CLUSTER:"+ "\n")
    fil.write("  - \"" + str(cluster) + "\"" +"\n")


for n in range(len(feature_list)):
    item = feature_list[n]
    path = feature_ref[n]
    command = 'mkdir ' + 'feature/' + str(item)
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,universal_newlines=True)
    command = 'scp '+ path +  ' feature/' + str(item) + "/"
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,universal_newlines=True)
    if args.bed is not None:
        if fbed[n] != "NA":
            bed_path = fbed[n]
            command = 'scp '+ bed_path +  ' feature/' + str(item) + "/"
            p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,universal_newlines=True)

if cluster == True:
    os.system('snakemake --cores ' + threads + ' all')
if cluster == False:
    os.system('snakemake all')
