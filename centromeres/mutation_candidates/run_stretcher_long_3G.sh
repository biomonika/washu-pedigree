#!/bin/bash -e
#SBATCH --nodes=1 
#SBATCH --job-name="run_stretcher_long.20250819"
#SBATCH --cpus-per-task=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=24G
#SBATCH --partition=long
#SBATCH --output=run_stretcher_long.20250819.%j.log

set -e
set -x

pwd; hostname; date

#ALIGN WITH STRETCHER
source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/groups/migalab/mcechova/conda/polishing

#submit with the grandparent file in the format chr*_extracted.fasta_grandparent.fa
gp=$1

# derive chromosome base
chrom="${gp%_grandparent.fa}"

mom="${chrom}_mother.fa"
gd="${chrom}_granddaughter.fa"

outdir="stretcher_long.${chrom%.*}"
mkdir -p "$outdir/generation1" "$outdir/generation2"

echo "[$chrom] generation1: grandparent vs mother"
srun -c 1 stretcher "$gp" "$mom" -outfile=$outdir/generation1/${outdir}.aln &

echo "[$chrom] generation2: mother vs granddaughter"
srun -c 1 stretcher "$mom" "$gd" -outfile=$outdir/generation2/${outdir}.aln &

wait
echo "[$chrom] done"

echo "Done."