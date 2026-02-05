#!/bin/bash

if [ -z ${SLURM_CPUS_PER_TASK+x} ]; then
    cpus=32
else
    cpus=$SLURM_CPUS_PER_TASK
fi

SCRIPT_DIR=$1
fastq=$2
ref=$3
out_dir=$4
group=$(basename $fastq | cut -d. -f1)

echo $group 

prefix=$out_dir/alignment/18S/$group
temp_prefix=$out_dir/alignment/18S/temp/$group

minimap2 \
    -t $cpus \
    -a \
    -x map-ont \
    -y \
    --MD \
    -Y \
    $ref \
    $fastq \
| samtools sort -@$cpus -O BAM --write-index -o $temp_prefix.nofilt.bam##idx##$temp_prefix.nofilt.bam.bai -
# | samtools sort -@$cpus -O BAM --write-index -o $temp_prefix.nofilter.bam##idx##$temp_prefix.nofilter.bam.bai -

samtools view -@$cpus -F 260 -h -O BAM --write-index -o $temp_prefix.bam##idx##$temp_prefix.bam.bai $temp_prefix.nofilt.bam

samtools view -@$cpus -h -O SAM -o $temp_prefix.sam $temp_prefix.bam 
samtools view -@$cpus $temp_prefix.bam > $temp_prefix.no_header.sam
samtools view -@$cpus -H $temp_prefix.bam > $temp_prefix.just_header.sam

paftools.js sam2paf $temp_prefix.sam > $temp_prefix.paf

python $SCRIPT_DIR/filter_both_18S.py $temp_prefix.no_header.sam $temp_prefix.paf sam 95 95 > $temp_prefix.filtered.sam
python $SCRIPT_DIR/filter_both_18S.py $temp_prefix.no_header.sam $temp_prefix.paf paf 95 95 > $temp_prefix.filtered.paf

awk '{print $1, $3, $4, $5}' $temp_prefix.filtered.paf > $temp_prefix.filtered.bed
bash $SCRIPT_DIR/get_chimeras_for_python.sh $temp_prefix.filtered.bed > $temp_prefix.chimeras.txt

if [[ $(cat $temp_prefix.chimeras.txt | wc -l) -eq 0 ]]; then
    cp $temp_prefix.filtered.paf $prefix.paf
    cp $temp_prefix.filtered.bed $prefix.bed
    cat $temp_prefix.just_header.sam $temp_prefix.filtered.sam | samtools sort -@$cpus -O BAM --write-index -o $prefix.bam##idx##$prefix.bam.bai -
else
    python $SCRIPT_DIR/remove_chimeras.py $temp_prefix.filtered.paf $temp_prefix.chimeras.txt > $prefix.paf
    python $SCRIPT_DIR/remove_chimeras.py $temp_prefix.filtered.bed $temp_prefix.chimeras.txt > $prefix.bed
    cat $temp_prefix.just_header.sam $temp_prefix.filtered.sam > $temp_prefix.temp.sam
    python $SCRIPT_DIR/remove_chimeras.py $temp_prefix.temp.sam $temp_prefix.chimeras.txt | samtools sort -@$cpus -O BAM --write-index -o $prefix.bam##idx##$prefix.bam.bai -
fi

num_alns=$(cat $prefix.paf | wc -l)
echo -e "final alignments (reads to ref) $num_alns"

if [[ $num_alns -gt 0 ]]; then
    samtools depth \
        -@$cpus \
        -a \
        $prefix.bam \
        | awk '{sum+=$3} END {print sum/NR}' > $prefix.18S_cov.txt
fi

rm -r $out_dir/alignment/18S/temp/$group.*
rm -r $out_dir/alignment/18S/$group.paf
rm -r $out_dir/alignment/18S/$group.bed
