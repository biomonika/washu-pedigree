# Analysis of rDNA G4 Stability and Sequence Conservation

This workflow analyzes the relationship between predicted G-quadruplex (G4) stability and sequence conservation across maternal and paternal rDNA assemblies. 

The workflow includes the following steps:
1. Download rDNA assemblies.
2. Annotate rDNA units using **andro**.
3. Multiple sequence alignment and consensus generation.
4. Re-aligning the sequences, and including the consensus.
5. Calculating the sequence conservation scores.
6. Processing G4Hunter predictions.
7. Comparing the G4 stability with the sequence conservation.
8. Visualizing the conservation across the G4 stability quartiles.

---

## 1. Download rDNA assemblies

The maternal and paternal rDNA sequences for PAN027 are downloaded from the assembly repository.

```
wget https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/rDNA/PAN027.rDNA.chr14.maternal.ref.fa
wget https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/rDNA/PAN027.rDNA.chr14.paternal.ref.fa
```

## 2. Annotate rDNA units using andro

Activate the environment containing andro and install the package.

`conda activate andro`
`pip install andro-0.3.tar.gz`

The rDNA assemblies are annotated using andro to identify rDNA repeat units.
```
andro PAN027.rDNA.chr14.maternal.ref.fa > PAN027.rDNA.chr14.maternal.ref.bed
andro PAN027.rDNA.chr14.paternal.ref.fa > PAN027.rDNA.chr14.paternal.ref.bed
```

Only annotated whole rDNA repeat units are retained for downstream analysis.

```
cat PAN027.rDNA.chr14.maternal.ref.bed | grep "RDNA" > PAN027.rDNA.chr14.maternal.RDNA.bed
cat PAN027.rDNA.chr14.paternal.ref.bed | grep "RDNA" > PAN027.rDNA.chr14.paternal.RDNA.bed
```

The corresponding sequences are extracted using bedtools getfasta.
```
bedtools getfasta \
    -fi PAN027.rDNA.chr14.maternal.ref.fa \
    -bed PAN027.rDNA.chr14.maternal.RDNA.bed \
    > PAN027.rDNA.chr14.maternal.RDNA.fa

bedtools getfasta \
    -fi PAN027.rDNA.chr14.paternal.ref.fa \
    -bed PAN027.rDNA.chr14.paternal.RDNA.bed \
    > PAN027.rDNA.chr14.paternal.RDNA.fa
```

## 3. Multiple sequence alignment and consensus generation

Install required alignment tools:

`conda install -c bioconda mafft emboss`

The extracted rDNA units are aligned independently for maternal and paternal haplotypes.


```
mafft --auto --maxiterate 100 PAN027.rDNA.chr14.paternal.RDNA.fa \
    > PAN027.rDNA.chr14.paternal.MSA.fa

mafft --auto --maxiterate 100 PAN027.rDNA.chr14.maternal.RDNA.fa \
    > PAN027.rDNA.chr14.maternal.MSA.fa
```

Consensus sequences are generated using EMBOSS cons.

Note: cons represents alignment gaps as N bases in the consensus sequence.

```
cons \
    -sequence PAN027.rDNA.chr14.maternal.MSA.fa \
    -outseq PAN027.rDNA.chr14.maternal.MSA.consensus.fa \
    -auto

cons \
    -sequence PAN027.rDNA.chr14.paternal.MSA.fa \
    -outseq PAN027.rDNA.chr14.paternal.MSA.consensus.fa \
    -auto
```

Rename consensus FASTA headers:

```
sed -i '' 's/^>EMBOSS_001/>maternal.MSA.consensus/g' \
    PAN027.rDNA.chr14.maternal.MSA.consensus.fa

sed -i '' 's/^>EMBOSS_001/>paternal.MSA.consensus/g' \
    PAN027.rDNA.chr14.paternal.MSA.consensus.fa
```

## 4. Re-aligning the sequences, and including the consensus

The consensus sequence is added back into the alignment and the MSA is recalculated.

```
cat PAN027.rDNA.chr14.maternal.MSA.consensus.fa \
    PAN027.rDNA.chr14.maternal.RDNA.fa \
    > PAN027.rDNA.chr14.maternal.RDNA+consensus.fa

cat PAN027.rDNA.chr14.paternal.MSA.consensus.fa \
    PAN027.rDNA.chr14.paternal.RDNA.fa \
    > PAN027.rDNA.chr14.paternal.RDNA+consensus.fa
```

Generate final alignments:

```
mafft --auto --maxiterate 100 PAN027.rDNA.chr14.paternal.RDNA+consensus.fa \
    > PAN027.rDNA.chr14.paternal.MSA.fa

mafft --auto --maxiterate 100 PAN027.rDNA.chr14.maternal.RDNA+consensus.fa \
    > PAN027.rDNA.chr14.maternal.MSA.fa
```

## 5. Calculating the sequence conservation scores.

Conservation scores are calculated across the final multiple sequence alignments. For each position in the alignment, the number of the basepairs matching the consensus is generated. 

```
python calculate_conservation.py \
    PAN027.rDNA.chr14.maternal.MSA.fa \
    PAN027.rDNA.chr14.maternal.MSA.conservation.bed \
    maternal.MSA.consensus

python calculate_conservation.py \
    PAN027.rDNA.chr14.paternal.MSA.fa \
    PAN027.rDNA.chr14.paternal.MSA.conservation.bed \
    paternal.MSA.consensus
```

## 6. Processing G4Hunter predictions.

G4Hunter predictions are exported from:

https://bioinformatics.ibp.cz/#/results/quadruplex

To run the analysis, upload the consensus rDNA sequences `PAN027.rDNA.chr14.maternal.MSA.consensus.fa` and `PAN027.rDNA.chr14.paternal.MSA.consensus.fa`.

Save the G4Hunter results as `PAN027.rDNA.chr14.maternal.MSA.consensus.csv` and `PAN027.rDNA.chr14.paternal.MSA.consensus.csv`.

Important: grouped.csv coordinates are 1-based, whereas BED coordinates are 0-based. We are interested in the grouped coordinates. The G4Hunter output is converted into BED format and assigned quartiles.

```
python reformat_G4Hunter_csv_to_bed_with_quartiles.py \
    PAN027.rDNA.chr14.maternal.MSA.consensus.csv \
    PAN027.rDNA.chr14.maternal.MSA.consensus.G4.bed \
    1.2 \
    maternal.MSA.consensus

python reformat_G4Hunter_csv_to_bed_with_quartiles.py \
    PAN027.rDNA.chr14.paternal.MSA.consensus.csv \
    PAN027.rDNA.chr14.paternal.MSA.consensus.G4.bed \
    1.2 \
    paternal.MSA.consensus
```

## 7. Comparing the G4 stability with the sequence conservation

Conservation scores are assigned to G4 regions using bedtools map.

The minimum conservation score overlapping each G4 region is reported.

*Maternal haplotype*
```
bedtools map \
    -a PAN027.rDNA.chr14.maternal.MSA.consensus.G4.bed \
    -b PAN027.rDNA.chr14.maternal.MSA.conservation.bed \
    -c 5 \
    -o min \
    > PAN027.rDNA.chr14.maternal.MSA.consensus.G4_vs_conservation.bed
```

*Paternal haplotype*
```
bedtools map \
    -a PAN027.rDNA.chr14.paternal.MSA.consensus.G4.bed \
    -b PAN027.rDNA.chr14.paternal.MSA.conservation.bed \
    -c 5 \
    -o min \
    > PAN027.rDNA.chr14.paternal.MSA.consensus.G4_vs_conservation.bed
```

## 8. Visualizing the conservation across the G4 stability quartiles

The final plots compare conservation scores across G4 score quartiles.
```
Rscript plot_conservation_boxplot.R \
    PAN027.rDNA.chr14.maternal.MSA.consensus.G4_vs_conservation.bed

Rscript plot_conservation_boxplot.R \
    PAN027.rDNA.chr14.paternal.MSA.consensus.G4_vs_conservation.bed
```
The resulting boxplots show that G-quadruplexes with higher stability were more likely to be located in variable rDNA regions (across units of the same haplotype).
