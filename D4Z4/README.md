## Resources for the annotation and plotting of the D4Z4 repeat unit 

Hailey Loucks
2025.12.02

### Annotation of the D4Z4 Repeat

HMM file - D4Z4rep.hmm
Built on CHM13 chr4:193538193-193541504
HMMer version 3.4

`nhmmer --cpu 12 --notextw --noali --tblout file.out  D4Z4rep.hmm genome.fa`

HMMer output files were converted to bed files using `./hmmertblout2bed.awk`

`gawk -v th=0.7 -f hmmertblout2bed.awk hmmer.out > hmmer.bed`

### Annotation of the pLAM repeat 

pLAM sequence obtained from Yeetong et. al. (PMID 37320968) 

Blast databases were constructed from the assemblies using  
`makeblastdb -in genome.fa -dbtype nucl -out genome` 
Blastn search performed with  
`blastn -db genomeDB -outfmt 6 -query pLAM.fa -num_threads 8 -out genome.out`

Functional poly-Adenylation (pLAM signal) identified manually by screening for presence of ATTAAA sequence (nonfunctional allele ATCAAA)