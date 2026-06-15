#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name="flag_3G_centromeres.20260608"
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=medium
#SBATCH --output=flag_3G_centromeres.20260608.%j.log

set -e
set -x 

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/.conda/envs/methylation

# analyze centromeres transmitted from the grandparents to the mother, that are also 
# passed to the granddaughter 
SAMPLES="PAN010 PAN011 PAN027 PAN028"

# prepare the files
cat transmissions.3G.txt | tail -n +2 | sed 's/\t/\n/g' > transmitted_centromeres.3G.txt

# download the needed annotations
for sample in $SAMPLES; do
    # download flagger CONSERVATIVE for centromeres
    wget -nc "https://public.gi.ucsc.edu/~mcechova/pedigree/annotations/v1.1/${sample}/flagger/flagger.${sample}.hifi.polished.v1.1.conservative.bed"
    wget -nc "https://public.gi.ucsc.edu/~mcechova/pedigree/annotations/v1.1/${sample}/flagger/flagger.${sample}.ONT.polished.v1.1.conservative.bed"
    #download nucflag
    wget -nc "https://public.gi.ucsc.edu/~mcechova/pedigree/annotations/v1.1/${sample}/nucflag.hifi.${sample}.bed"
    #download active_hor
    wget -nc "https://public.gi.ucsc.edu/~mcechova/pedigree/annotations/v1.1/${sample}/${sample}.polished.cenSat.active_hor.bed"
done

touch bed.md5sum.txt
for a in *.bed; do
    md5sum "$a" >> bed.md5sum.txt
done

############################################
# SUM FLAGGER HIFI ANNOTATIONS INSIDE CENTROMERES
############################################

for sample in $SAMPLES; do
    centromere_bed="${sample}.polished.cenSat.active_hor.bed"
    flagger_bed="flagger.${sample}.hifi.polished.v1.1.conservative.bed"

    if [[ ! -f "$centromere_bed" || ! -f "$flagger_bed" ]]; then
        echo "Skipping ${sample}: missing $centromere_bed or $flagger_bed"
        continue
    fi

    bedtools intersect \
        -a "$centromere_bed" \
        -b "$flagger_bed" \
        -wo \
    | cut -f1,13,19 \
    | sort \
    | awk '{ key = $1 "\t" $2; sum[key] += $3 } END { for (k in sum) print k "\t" sum[k] }' \
    | sort \
    > ${sample}.flagger_hifi.problems_in_centromeres.txt

    grep -F -f transmitted_centromeres.3G.txt \
        ${sample}.flagger_hifi.problems_in_centromeres.txt \
        | sort -V \
        > transmitted_flagger_hifi_issues_in_${sample}_centromeres.txt
done

############################################
# SUM FLAGGER ONT ANNOTATIONS INSIDE CENTROMERES
############################################

for sample in $SAMPLES; do
    centromere_bed="${sample}.polished.cenSat.active_hor.bed"
    flagger_bed="flagger.${sample}.ONT.polished.v1.1.conservative.bed"

    if [[ ! -f "$centromere_bed" || ! -f "$flagger_bed" ]]; then
        echo "Skipping ${sample}: missing $centromere_bed or $flagger_bed"
        continue
    fi

    bedtools intersect \
        -a "$centromere_bed" \
        -b "$flagger_bed" \
        -wo \
    | cut -f1,13,19 \
    | sort \
    | awk '{ key = $1 "\t" $2; sum[key] += $3 } END { for (k in sum) print k "\t" sum[k] }' \
    | sort \
    > ${sample}.flagger_ONT.problems_in_centromeres.txt

    grep -F -f transmitted_centromeres.3G.txt \
        ${sample}.flagger_ONT.problems_in_centromeres.txt \
        | sort -V \
        > transmitted_flagger_ONT_issues_in_${sample}_centromeres.txt
done

############################################
# SUM NUCFLAG ANNOTATIONS INSIDE CENTROMERES
############################################

for sample in $SAMPLES; do
    centromere_bed="${sample}.polished.cenSat.active_hor.bed"
    nucflag_bed="nucflag.hifi.${sample}.bed"

    if [[ ! -f "$centromere_bed" || ! -f "$nucflag_bed" ]]; then
        echo "Skipping ${sample}: missing $centromere_bed or $nucflag_bed"
        continue
    fi

    bedtools intersect \
        -a "$centromere_bed" \
        -b "$nucflag_bed" \
        -wo \
    | cut -f1,13-14 \
    | sort \
    | awk '{ key = $1 "\t" $2; sum[key] += $3 } END { for (k in sum) print k "\t" sum[k] }' \
    | sort \
    > ${sample}.nucflag.problems_in_centromeres.txt

    grep -F -f transmitted_centromeres.3G.txt \
        ${sample}.nucflag.problems_in_centromeres.txt \
        | sort -V \
        > transmitted_nucflag_hifi_issues_in_${sample}_centromeres.txt
done

echo "Done."

# manually write these into our spreadsheet