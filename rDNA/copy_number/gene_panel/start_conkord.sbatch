#!/usr/bin/bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=65G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/%A.out
#SBATCH --error=logs/%A.err

ml python
ml jellyfish
ml bedtools
ml samtools

python conkord.py -k 31 -bed gene_fasta_files/C13orf46.bed gene_fasta_files/CCDC85C.bed gene_fasta_files/LINC00322.bed gene_fasta_files/LINC02205.bed gene_fasta_files/PGF.bed gene_fasta_files/REM2.bed gene_fasta_files/CACNA1I.bed gene_fasta_files/HIC2.bed gene_fasta_files/INAFM2.bed gene_fasta_files/LINC01656.bed gene_fasta_files/PLA2G4F.bed gene_fasta_files/SNX33.bed -f gene_fasta_files/C13orf46.fa gene_fasta_files/CCDC85C.fa gene_fasta_files/LINC00322.fa gene_fasta_files/LINC02205.fa gene_fasta_files/PGF.fa gene_fasta_files/REM2.fa gene_fasta_files/CACNA1I.fa gene_fasta_files/HIC2.fa gene_fasta_files/INAFM2.fa gene_fasta_files/LINC01656.fa gene_fasta_files/PLA2G4F.fa gene_fasta_files/SNX33.fa -r sequencing_data/ -t 15 --cluster -g /n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg23/gene_normalization_panel/chm13_renamed.fa --gzip
