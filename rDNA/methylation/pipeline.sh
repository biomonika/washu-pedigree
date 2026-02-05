#!/bin/bash

ml load modkit

SCRIPT_DIR=$1
bam=$2
chr_list_f=$3
ref=$4
out_dir=$5
genome_bed=$6
eighteen_S_fa=$7
fortyfive_S_bed=$8
TR_fa=$9

map_jids=""
meth_jids=""
TR_jids=""
eighteen_S_jids=""

mkdir -p $out_dir $out_dir/logs $out_dir/fastqs
mkdir -p $out_dir/alignment $out_dir/alignment/logs $out_dir/alignment/temp $out_dir/alignment/read_breakdown 
mkdir -p $out_dir/alignment/18S $out_dir/alignment/18S/logs $out_dir/alignment/18S/temp
mkdir -p $out_dir/get_methylation $out_dir/get_methylation/modkit_beds $out_dir/get_methylation/logs $out_dir/get_methylation/read_breakdown
mkdir -p $out_dir/final_analysis

for fastq in $out_dir/fastqs/*.fq.gz; do
	echo $fastq
    group=$(basename $fastq | cut -d. -f1)
    echo "calling on $group"
	# if [ ! -s $out_dir/alignment/logs/$group.success.txt ]; then
    	map_jid="$(sbatch \
			--time=1-00:00:00 \
			--partition=norm \
			--mem=16g \
			--cpus-per-task=16 \
			-o $out_dir/alignment/logs/$group.out \
			$SCRIPT_DIR/map_to_ref_and_split_per_group.sh \
				$SCRIPT_DIR \
				$fastq \
				$ref \
				$out_dir)"
		map_jids+=$map_jid

		echo "$group map job id : $map_jid"
	# fi

	eighteen_S_jid="$(sbatch \
				--time=01:00:00 \
				--mem=16g \
				--partition=norm,quick \
				--cpus-per-task=8 \
				-o $out_dir/alignment/18S/logs/$group.out \
				$SCRIPT_DIR/get_18S_cov.sh \
					$SCRIPT_DIR \
					$fastq \
					$eighteen_S_fa \
					$out_dir)"

	echo "$group 18S job id : $eighteen_S_jid"
	eighteen_S_jids+=":$eighteen_S_jid"

	if [[ ${#map_jid} -gt 0 ]]; then
		meth_jid="$(sbatch \
						--time=01:00:00 \
						--partition=norm,quick \
						--mem=4g \
						-o $out_dir/get_methylation/logs/$group.whole_group.out \
						--dependency=afterok:$map_jid \
						$SCRIPT_DIR/get_methylation_whole_group.sh \
							$ref \
							$group \
							$out_dir \
							$fortyfive_S_bed)"

		echo "$group methylation job id : $meth_jid"
		meth_jids+=":$meth_jid"
	else
		meth_jid="$(sbatch \
						--time=01:00:00 \
						--partition=norm,quick \
						--mem=4g \
						-o $out_dir/get_methylation/logs/$group.whole_group.out \
						$SCRIPT_DIR/get_methylation_whole_group.sh \
							$ref \
							$group \
							$out_dir \
							$fortyfive_S_bed)"
		echo "$group methylation job id : $meth_jid"
		meth_jids+=":$meth_jid"
	fi

	if [[ $# -gt 7 ]]; then # get methylation breakdown of the subregion
		if [[ ${#map_jid} -gt 0 ]]; then
			meth_jid="$(sbatch \
						--time=08:00:00 \
						--mem=8g \
						--cpus-per-task=16 \
						-o $out_dir/get_methylation/logs/$group.breakdown.out \
						--dependency=afterok:$map_jid \
						$SCRIPT_DIR/get_methylation_breakdown.sh \
							$ref \
							$group \
							$out_dir \
							$fortyfive_S_bed)"

			echo "$group methylation breakdown job id : $meth_jid"
			meth_jids+=":$meth_jid"
		else
			meth_jid="$(sbatch \
						--time=08:00:00 \
						--mem=8g \
						--cpus-per-task=16 \
						-o $out_dir/get_methylation/logs/$group.breakdown.out \
						$SCRIPT_DIR/get_methylation_breakdown.sh \
							$ref \
							$group \
							$out_dir \
							$fortyfive_S_bed)"
			echo "$group methylation breakdown job id : $meth_jid"
			meth_jids+=":$meth_jid"
		fi
	fi

	if [[ $# -eq 9 ]]; then # get TR copy num
		if [[ ${#map_jid} -gt 0 ]]; then
			TR_jid="$(sbatch \
							--time=04:00:00 \
							--partition=norm,quick \
							--mem=8g \
							-o $out_dir/get_TR/logs/$group.out \
							--dependency=afterok:$map_jid \
							$SCRIPT_DIR/get_TR_breakdown.sh \
								$SCRIPT_DIR \
								$group \
								$out_dir \
								$TR_fa)"
			echo "$group TR job id : $TR_jid"
			TR_jids+=":$TR_jid"
		else
			TR_jid="$(sbatch \
							--time=04:00:00 \
							--partition=norm,quick \
							--mem=8g \
							-o $out_dir/get_TR/logs/$group.out \
							$SCRIPT_DIR/get_TR_breakdown.sh \
								$SCRIPT_DIR \
								$group \
								$out_dir \
								$TR_fa)"
			echo "$group TR job id : $TR_jid"
			TR_jids+=":$TR_jid"
		fi
	fi
done

all_jids=""
all_jids+=$map_S_jids
all_jids+=$eighteen_S_jids
all_jids+=$meth_jids
all_jids+=$TR_jids

if [[ ${all_jids:0:1} == ":" ]]; then
	all_jids="${all_jids:1}"
fi

all_jids=$(echo $all_jids | awk '{sub(/:$/, "", $0); print}')

echo "all dependency job ids - $all_jids"

if [[ ${#all_jids} -gt 0 ]]; then
	if [[ $# -eq 9 ]]; then
		sbatch --time=30:00 --mem=4g --partition=norm,quick -o $out_dir/logs/summary.txt --dependency=afterok:$all_jids $SCRIPT_DIR/summary_and_analysis.sh $SCRIPT_DIR $out_dir 1
	else
	if [[ -z "$fortyfive_S_bed" ]]; then
		sbatch --time=30:00 --mem=4g --partition=norm,quick -o $out_dir/logs/summary.txt --dependency=afterok:$all_jids $SCRIPT_DIR/summary_and_analysis.sh $SCRIPT_DIR $out_dir
	else
		sbatch --time=30:00 --mem=4g --partition=norm,quick -o $out_dir/logs/summary.txt --dependency=afterok:$all_jids $SCRIPT_DIR/summary_and_analysis.sh $SCRIPT_DIR $out_dir 1
	fi
else
	if [[ $# -gt 8 ]]; then
		sbatch --time=30:00 --mem=4g --partition=norm,quick -o $out_dir/logs/summary.txt $SCRIPT_DIR/summary_and_analysis.sh $SCRIPT_DIR $out_dir 1
	else
	if [[ -z "$fortyfive_S_bed" ]]; then
		sbatch --time=30:00 --mem=4g --partition=norm,quick -o $out_dir/logs/summary.txt $SCRIPT_DIR/summary_and_analysis.sh $SCRIPT_DIR $out_dir
	else
		sbatch --time=30:00 --mem=4g --partition=norm,quick -o $out_dir/logs/summary.txt $SCRIPT_DIR/summary_and_analysis.sh $SCRIPT_DIR $out_dir 1
	fi
fi

sbatch --time=30:00 --mem=4g --partition=norm,quick --dependency=afterok:$all_jids -o $out_dir/logs/summary.txt $SCRIPT_DIR/summarize_group_meth.sh $SCRIPT_DIR $out_dir
