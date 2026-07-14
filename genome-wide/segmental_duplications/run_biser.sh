#!/bin/bash
#SBATCH --job-name=run_biser.20241122
#SBATCH --partition=long
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --output=run_biser.20241122.%j.log
#SBATCH --time=12:00:00

set -e
set -x

pwd; hostname; date

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/groups/migalab/mcechova/conda/biser

assembly=$1
assembly_name="${assembly%.*}"

biser -o biser.${assembly_name} --keep-contigs -t 2 --gc-heap 10G ${assembly}
echo "Done."
date

