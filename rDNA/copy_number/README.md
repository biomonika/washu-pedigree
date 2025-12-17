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

# Washu Processing

1.  The pipeline CONKORD was used to estimate rDNA copy numbers for the WashU pedigree. It was launched on a cluster with the script "start_conkord.sbatch". The invocation was:
python conkord.py --no_uniq -k 31 -bed rdnaModel.bed -f rdna/KY962518_18s.fa -r sequencing_data/ -t 15 --cluster -g chm13v2.0_singleline_unmasked.fa --gzip
2.  Base rDNA copy number estimates were further normalized to a panel of single-copy genes with similar G/C content to the 18S subunit. These genes were also estimated with CONKORD (altered to produce non-integer outputs). The expected diploid copy number for each gene is 2, so for each sample the panel average copy number was divided by 2, then that factor was multiplied by the sample rDNA copy number to adjust for expected over or under estimation. Genes used were:
C13orf46, CACNA1I, CCDC85C, HIC2, IGHV169D, INAFM2, LINC00322, LINC01656, LINC02205, MIR6882, PGF, PLA2G4F, REM2, SNX33, TRAJ18
