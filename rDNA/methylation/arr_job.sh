#!/bin/bash

ml load bedtools minimap2 samtools seqtk
mkdir -p logs

SCRIPT_DIR=/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/absolute_path_rdna_analysis_scripts

i=$SLURM_ARRAY_TASK_ID

if [[ -z $i ]]; then
  echo "No slurm array task ID nor line_num given. Exit"
  exit -1
fi

sample=`sed -n ${i}p samples.txt | awk '{print $1}'`

bash $SCRIPT_DIR/pipeline.sh \
					$SCRIPT_DIR \
					fake.bam \
					fake_chrs.txt \
					/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/KY962518-ROT.fa \
					out_45S/$sample \
					fake.genome.bed \
					/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/18S_on_KY-ROT.fa \
					/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/45S_on_KY-ROT.bed \
					/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/TR.fa

# mkdir -p out_promoter/$sample out_promoter/$sample/fastqs
# cp out_45S/$sample/fastqs/* out_promoter/$sample/fastqs/
# bash $SCRIPT_DIR/pipeline.sh \
# 					$SCRIPT_DIR \
# 					fake.bam \
# 					fake_chrs.txt \
# 					/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/KY962518-ROT.fa \
# 					out_promoter/$sample \
# 					fake.genome.bed \
# 					/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/18S_on_KY-ROT.fa \
# 					/data/Phillippy2/projects/chm13_rdna_methylation_reanalysis/refs/rdna_unit/main_promoter.bed
