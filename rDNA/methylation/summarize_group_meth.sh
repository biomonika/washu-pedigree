#!/bin/bash

SCRIPT_DIR=$1
out_dir=$2

for bam in $out_dir/alignment/*.bam; do
	group=$(basename $bam | cut -d. -f1)
	avg_meth_pct=$(awk '{sum+=($11/100)} END {print sum/NR}' $out_dir/get_methylation/modkit_beds/$group.whole_region.bed)
	echo -e "$group\t$avg_meth_pct" >> $out_dir/final_analysis/group_summary.txt
done