#!/bin/bash

ml load modkit 

if [ -z ${SLURM_CPUS_PER_TASK+x} ]; then
    cpus=32
else
    cpus=$SLURM_CPUS_PER_TASK
fi

get_meth_breakdown_read()
{
	ref=$1
	group=$2
	out_dir=$3
	readname=$4
	subregion_bed=$5

	path=$out_dir/get_methylation/read_breakdown/$group/$readname
	mkdir -p $path $path/modkit_beds $path/logs

	for aln in $out_dir/alignment/read_breakdown/$group/$readname/*.bam; do
		num=$(basename $aln | cut -d. -f1)
		if [ -z "$subregion_bed" ]; then
			echo "$readname $num no subregion"
			modkit pileup \
				--threads $cpus \
				--ref $ref \
				--cpg \
				--combine-strands \
				--ignore h \
				--log $path/logs/$num.log \
				$aln \
				$path/modkit_beds/$num.bed
		else
			subregion=$(awk '{printf "%s:%d-%d", $1, $2, $3}' $subregion_bed)
			echo "$readname $num subregion $subregion"
			modkit pileup \
				--threads $cpus \
				--ref $ref \
				--region $subregion \
				--cpg \
				--combine-strands \
				--ignore h \
				--log $path/logs/$num.log \
				$aln \
				$path/modkit_beds/$num.bed
		fi
	done
}

MAX_JOBS=$cpus

ref=$1
group=$2
out_dir=$3
subregion_bed=$4

echo "
ref=$1
group=$2
out_dir=$3
subregion_bed=$4
"

echo $group

mkdir -p $out_dir/get_methylation/read_breakdown/$group

if [[ $(ls $out_dir/alignment/read_breakdown/$group/* | wc -l) -gt 0 ]]; then
	for fp in $out_dir/alignment/read_breakdown/$group/*; do
		readname=$(basename $fp)
		while [ `jobs -p | wc -l` -ge ${MAX_JOBS} ]
		do
			sleep 5
		done
		get_meth_breakdown_read $ref $group $out_dir $readname $subregion_bed &
	done
	wait `jobs -p`
fi
