# CONKORD Usage
A snakemake pipeline for counting copy number of genomic features. The python file conkord.py creates a config.yml file for snakemake then runs a "Snakemake all" command.

Example Usage:

python conkord.py --no_uniq -k 31 -bed feature_coordinates.bed -f feature.fa -r sequencing_data/ -t 15 --cluster -g genome.fa --gzip  

Parameter Description:  
--no_uniq:  Do not use only unique kmers for your feature. This is off by default but is useful for features like rDNA where it is hard to identify all of the genomic locations.  
-k:  kmer size (default 31)  
-f:  The fasta for your feature of interest. It must be a SINGLE LINE fasta (as in one line per fasta entry).  
-bed  A bed file with coordinates for your feature of interest in the reference genome.  
-r:  Directory with only the sequencing reads in fastq or fastq.gz format. Reads must be paired end and follow a name scheme ending with _1.fastq and _2.fastq or _1.fastq.gz and _2.fastq.gz (not R1 or R2).  
-g:  The reference genome. It must be a SINGLE LINE fasta (as in one line per fasta entry).  
-gzip:  If the reads are gzipped you must pass this flag  
-t:  Threads  
-w_size:  The window size for finding G/C matched windows to normalize to (default 2000)  
--cluster:  If you are using a cluster with slurm this will cause conkord.py to submit an sbatch command for snakemake with the number of cores provided  
