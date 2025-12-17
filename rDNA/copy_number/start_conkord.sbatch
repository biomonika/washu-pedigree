#!/usr/bin/bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=150G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/%A.out
#SBATCH --error=logs/%A.err

ml python
ml jellyfish
ml bedtools
ml samtools

python conkord.py --no_uniq -k 31 -bed rdnaModel.bed -f rdna/KY962518_18s.fa -r sequencing_data/ -t 15 --cluster -g /n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/chm13_assembly/v2.0/chm13v2.0_singleline_unmasked.fa --gzip
