#!/bin/bash

# Edit CHM13v2 genetic map files to format expected by flare

MAP=genetic_map_chm13
mkdir -p $MAP/old_maps

for i in {1..22} X
do
	file=$MAP/chr${i}_no_mask.txt
	# Replace comma with tabs
	# Reorder columns (chr, SNP ID, cumulative cM, base position)
	awk 'BEGIN {FS=","; OFS="\t"} NR>1 {print $1, ".", $6, $2}' $file \
		> $MAP/chr${i}.CHM13v2.map

	mv $file $MAP/old_maps
done