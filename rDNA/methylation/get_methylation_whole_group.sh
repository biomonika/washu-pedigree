#!/bin/bash

ml load modkit 

if [ -z ${SLURM_CPUS_PER_TASK+x} ]; then
    cpus=32
else
    cpus=$SLURM_CPUS_PER_TASK
fi

ref=$1
group=$2
out_dir=$3
subregion_bed=$4
echo $group

if [[ $(samtools view -@$cpus -c $out_dir/alignment/$group.bam) -eq 0 ]]; then
	echo "empty bam"
	exit 0
fi

if [[ -z "$subregion_bed" ]]; then
	modkit pileup \
		--threads $cpus \
		--ref $ref \
		--cpg \
		--combine-strands \
		--ignore h \
		--log $out_dir/get_methylation/logs/$group.whole_region.log \
		$out_dir/alignment/$group.bam \
		$out_dir/get_methylation/modkit_beds/$group.whole_region.bed
else
	modkit pileup \
		--threads $cpus \
		--ref $ref \
		--cpg \
		--combine-strands \
		--ignore h \
		--log $out_dir/get_methylation/logs/$group.whole_region.log \
		$out_dir/alignment/$group.bam \
		$out_dir/get_methylation/modkit_beds/$group.whole_region.bed

	subregion=$(awk '{printf "%s:%d-%d", $1, $2, $3}' $subregion_bed)
	modkit pileup \
		--threads $cpus \
		--ref $ref \
		--region $subregion \
		--cpg \
		--combine-strands \
		--ignore h \
		--log $out_dir/get_methylation/logs/$group.subregion.log \
		$out_dir/alignment/$group.bam \
		$out_dir/get_methylation/modkit_beds/$group.subregion.bed
fi
