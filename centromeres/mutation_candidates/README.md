# Identify candidate mutations in centromeres

The goal of this section is to identify the candidate three-generational and two-generational candidate mutations in centromeres. 

## Prepare the centromeric fasta files
`cat PAN*.fasta >pedigree.fasta`
`cat PAN*active_hor.bed >pedigree.active_hor.bed`

Which chromosomes have maternal and paternal transmissions?

`cat transmissions.txt | sed 's/\t/\n/g' | grep "PAN027" | grep "maternal" | cut -d'.' -f2 | sort >maternal.transmissions.txt`
`cat transmissions.txt | sed 's/\t/\n/g' | grep "PAN027" | grep "paternal" | cut -d'.' -f2 | sort >paternal.transmissions.txt`

Extract the centromeres into individual files, while merging nearby bed annotations to account for gaps

`sbatch extract_centromeres_3G.sh`

## 3 GENERATIONS

### Align active_hor centromeric regions using minimap2/paftools and stretcher

```
for a in chr*_grandparent.fa; do 
    echo $a; 
    sbatch run_stretcher_long_3G.sh $a; 
done;
```

Run chr3 individually since it has two active arrays, skip the problematic chr4.

`sbatch run_stretcher_long_3G.sh chr3_1_grandparent.fa`

`sbatch run_stretcher_long_3G.sh chr3_2_grandparent.fa`

Align one centromere from each generation

`for a in chr*.fasta; do echo $a; sbatch align_centromeres_3G.sh $a; done;`

Intersect to find out which variants from VCF file overlap active_hor (if using broad centromere definition from the CenSat annotation track)

`sbatch intersect_VCF_files_with_active_hor.sh`

Loop through all .vcf files and index (uncompressed)
```
for vcf in *active_hor.vcf; do
    echo "Compressing: $vcf"
    bgzip -f "$vcf"  # overwrite with .vcf.gz

    gzipped_vcf="${vcf}.gz"
    echo "Indexing: $gzipped_vcf"
    tabix -p vcf "$gzipped_vcf"
done
```

### Identify candidate de-novo mutations

List of chromosomes to process
```
chromosomes=(chr{1..22} chrX)
for chr in "${chromosomes[@]}"; do
	echo $chr
    vcf1=(${chr}_extracted.fasta_generations1_shifted_*aternal_active_hor.vcf.gz)
    vcf2=(${chr}_extracted.fasta_generations2_shifted_*aternal_active_hor.vcf.gz)

    echo $vcf1
    echo $vcf2

    output="${chr}.candidate_mutation.vcf"
    echo "Comparing: $vcf1 vs $vcf2"
    bcftools isec -C -w1 -O v -o "$output" "$vcf1" "$vcf2"
done
```
Plot those candidate mutations

`for a in *candidate_mutation.vcf; do echo $a; python plot_mutations_from_vcf.py --vcf ${a}; done;`

## 2 GENERATIONS

Run stretcher alignments for 2G

`for a in *grandparent.fa; do echo $a; sbatch run_stretcher_long.2G.sh $a; done;`

```
###########################################################
# find mutations in all chromosomes using minimap2+paftools
###########################################################

for chr in {1..22} X; do
    echo "Processing chromosome: chr${chr}"
    sbatch align_2G_centromeres.sh "chr"${chr};
done
```

### Parse the results

```
##############
#detect indels
##############

for a in *shifted.vcf; do
    awk -F'\t' '
        !/^#/ && length($4) != length($5) {
            $4 = (length($4) > 2 ? "LEN(" length($4)-1 ")" : $4)
            $5 = (length($5) > 2 ? "LEN(" length($5)-1 ")" : $5)
            print $1, $2, $4, $5
        }
    ' OFS='\t' "$a"
done | sort -V | grep -v "chr4"
```
```
##############
#detect SNVs
##############

for a in *shifted.vcf; do
    awk -F'\t' '
        !/^#/ && length($4) == 1 && length($5) == 1 {
            print $1, $2, $4, $5
        }
    ' OFS='\t' "$a"
done | sort -V | grep -v "chr4"
```
