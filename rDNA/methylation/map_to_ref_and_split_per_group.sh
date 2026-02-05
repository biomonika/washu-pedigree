#!/bin/bash

# set -m

if [ -z ${SLURM_CPUS_PER_TASK+x} ]; then
    cpus=32
else
    cpus=$SLURM_CPUS_PER_TASK
fi


split_on_char()
{
    SCRIPT_DIR=$1
    breakdown_dir=$2
    prefix=$3
    temp_prefix=$4
    fastq=$5
    first_char=$6

    echo "call split_on_char $SCRIPT_DIR $breakdown_dir $prefix $temp_prefix $fastq $first_char"
    # echo "SCRIPT_DIR=$1
    # breakdown_dir=$2
    # prefix=$3
    # temp_prefix=$4
    # fastq=$5
    # first_char=$6"

    while read line; do
    readname=$(echo $line | awk '{print $1}')
    # if [ -f $breakdown_dir/$readname/$readname.fa ]; then
    #     continue
    # fi
    flag=$(echo $line | awk '{print $2}')
    echo $readname
    mkdir -p $breakdown_dir/$readname
    echo $readname > $breakdown_dir/$readname/$readname.lst
    if [[ $flag -eq 16 ]]; then
            seqtk subseq $fastq $breakdown_dir/$readname/$readname.lst | seqtk seq -r | bgzip -@$cpus > $breakdown_dir/$readname/$readname.fq.gz
        echo "True" > $breakdown_dir/$readname/$readname.reverse
    else
            seqtk subseq $fastq $breakdown_dir/$readname/$readname.lst | bgzip -@$cpus > $breakdown_dir/$readname/$readname.fq.gz
        echo "False" > $breakdown_dir/$readname/$readname.reverse
    fi
    seqtk seq -a $breakdown_dir/$readname/$readname.fq.gz > $breakdown_dir/$readname/$readname.fa
    samtools faidx $breakdown_dir/$readname/$readname.fa
    done <  <(samtools view -@$cpus $prefix.bam | grep "^$first_char")

    samtools view -@$cpus $prefix.bam | grep "^$first_char" | python $SCRIPT_DIR/split.py $breakdown_dir $temp_prefix.just_header.sam

    if [[ $(ls $breakdown_dir/$first_char*/*.sam | wc -l) -gt 0 ]]; then
    for sam in $breakdown_dir/$first_char*/*.sam; do
        sam_path=$(echo $sam | cut -d. -f 1)
        readname=$(echo $sam | cut -d/ -f 5)
        readlen=$(awk '{print $2}' $breakdown_dir/$readname/$readname.fa.fai)
        samtools view -@$cpus -h -O BAM --write-index -o $sam_path.bam##idx##$sam_path.bam.bai $sam
        if [[ $(cat $breakdown_dir/$readname/$readname.reverse) == "True" ]]; then
            echo "$readname is reversed"
            paftools.js sam2paf $sam_path.sam | awk -v OFS='\t' -v len=$readlen '{print $1,len-$4,len-$3}' > $sam_path.bed
        else
            echo "$readname is already forward"
            paftools.js sam2paf $sam_path.sam | awk -v OFS='\t' '{print $1,$3,$4}' > $sam_path.bed
        fi
        bedtools getfasta -fi $breakdown_dir/$readname/$readname.fa -bed $sam_path.bed -fo $sam_path.fa
    done
    fi
}


echo $cpus 

SCRIPT_DIR=$1
fastq=$2
ref=$3
out_dir=$4
group=$(basename $fastq | cut -d. -f1)

echo "
SCRIPT_DIR=$1
fastq=$2
ref=$3
out_dir=$4"
echo "-----"
echo $group 

prefix=$out_dir/alignment/$group
temp_prefix=$out_dir/alignment/temp/$group

minimap2 \
    -t $cpus \
    -a \
    -x map-ont \
    -y \
    --MD \
    -Y \
    $ref \
    $fastq \
    | samtools view -@$cpus -F 260 -h -O SAM -o $temp_prefix.sam

samtools view -@$cpus $temp_prefix.sam > $temp_prefix.no_header.sam
samtools view -@$cpus -H $temp_prefix.sam > $temp_prefix.just_header.sam

paftools.js sam2paf $temp_prefix.sam > $temp_prefix.paf

python $SCRIPT_DIR/filter_both.py $temp_prefix.no_header.sam $temp_prefix.paf sam > $temp_prefix.filtered.sam
python $SCRIPT_DIR/filter_both.py $temp_prefix.no_header.sam $temp_prefix.paf paf > $temp_prefix.filtered.paf

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
    breakdown_dir=$out_dir/alignment/read_breakdown/$group
    mkdir -p $breakdown_dir
    while read first_char; do
        echo $first_char
        split_on_char $SCRIPT_DIR $breakdown_dir $prefix $temp_prefix $fastq $first_char &
    done < <(samtools view -@$cpus $prefix.bam | awk '{print substr($1,1,1)}' | sort -u) 
    wait `jobs -p`
fi

rm -r $out_dir/alignment/temp/$group.*
rm -r $out_dir/alignment/$group.paf
rm -r $out_dir/alignment/$group.bed

touch $out_dir/alignment/logs/$group.success.txt
