#!/bin/bash
#SBATCH --job-name=process_biser.20241122
#SBATCH --partition=long
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --output=process_biser.20241122.%j.log
#SBATCH --time=12:00:00

set -e
set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/groups/migalab/mcechova/conda/biser

#the default output from biser is .bedpe format
for f in *softmasked; do
    mv "$f" "${f}.bedpe"
done

#reformat bedpe to bed format
for a in biser*.bedpe; do
  base="${a%.bedpe}"
  awk '{
    name = ($7 ? $7 : ".")
    score = ($8 ? $8 : "0")
    strand1 = ($9 ? $9 : ".")
    strand2 = ($10 ? $10 : ".")
    color = "128,128,128"
    print $1, $2, $3, name, score, strand1, $2, $3, color
    print $4, $5, $6, name, score, strand2, $5, $6, color
  }' OFS="\t" "$a" > "${base}.bed"
  echo "✓ Converted $a → ${base}.bed"
done

date

