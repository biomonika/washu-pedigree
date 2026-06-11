# This script prepares the annotations for chromosome-level plotting 
# of meiotic recombination breakpoints, ancestry, 
# and repeats/satellites (cenSat annotation track)

set -e
set -x

#################################################################
# format the ancestry bed files
#################################################################

#skip the first line because it includes the header
tail -n +2 PAN027.ancestry.maternal.assemblyCoords.bed | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' | bedtools sort >PAN027.ancestry.maternal.bed
tail -n +2 PAN027.ancestry.paternal.assemblyCoords.bed | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' | bedtools sort >PAN027.ancestry.paternal.bed

#################################################################
# MATERNAL lineage, prepare the data for the meiotic recombination breakpoints
#################################################################

cat PAN010_sharedcontigs_to_PAN027-hap1.bed | awk '{print $1"\t"$3"\t"$4}' | bedtools merge | awk '{print $0"\t#9EE2BF"}' >maternal.hap1.breakpoints.bed
cat PAN010_sharedcontigs_to_PAN027-hap2.bed | awk '{print $1"\t"$3"\t"$4}' | bedtools merge | awk '{print $0"\t#6495EC"}' >maternal.hap2.breakpoints.bed

#remove potential overlaps between hap1 and hap2
bedtools subtract -a maternal.hap1.breakpoints.bed -b maternal.hap2.breakpoints.bed >tmp.maternal.hap1.breakpoints.bed
bedtools subtract -a maternal.hap2.breakpoints.bed -b maternal.hap1.breakpoints.bed >tmp.maternal.hap2.breakpoints.bed

cat tmp.maternal.hap1.breakpoints.bed tmp.maternal.hap2.breakpoints.bed | bedtools sort >maternal.breakpoints.bed
rm tmp.maternal.hap1.breakpoints.bed tmp.maternal.hap2.breakpoints.bed

#################################################################
# PATERNAL lineage, prepare the data for the meiotic recombination breakpoints
#################################################################

cat PAN011_sharedcontigs_to_PAN027-hap1.bed | awk '{print $1"\t"$3"\t"$4}' | bedtools merge | awk '{print $0"\t#9EE2BF"}' >paternal.hap1.breakpoints.bed
cat PAN011_sharedcontigs_to_PAN027-hap2.bed | awk '{print $1"\t"$3"\t"$4}' | bedtools merge | awk '{print $0"\t#6495EC"}' >paternal.hap2.breakpoints.bed

#remove potential overlaps between hap1 and hap2
bedtools subtract -a paternal.hap1.breakpoints.bed -b paternal.hap2.breakpoints.bed >tmp.paternal.hap1.breakpoints.bed
bedtools subtract -a paternal.hap2.breakpoints.bed -b paternal.hap1.breakpoints.bed >tmp.paternal.hap2.breakpoints.bed

cat tmp.paternal.hap1.breakpoints.bed tmp.paternal.hap2.breakpoints.bed | bedtools sort >paternal.breakpoints.bed
rm tmp.paternal.hap1.breakpoints.bed tmp.paternal.hap2.breakpoints.bed

#################################################################
# MERGE THE CONSECUTIVE ANNOTATIONS
#################################################################

for hap in maternal paternal; do
    awk '
    BEGIN { OFS="\t" }
    {
      if (NR == 1 || $4 != prev_color || $1 != prev_chr) {
        if (NR > 1) print prev_chr, start, prev_end, prev_color;
        start = $2;
      }
      prev_chr = $1;
      prev_end = $3;
      prev_color = $4;
    }
    END {
      print prev_chr, start, prev_end, prev_color;
    }' "${hap}.breakpoints.bed" > "${hap}.breakpoints.merged.bed"
done



