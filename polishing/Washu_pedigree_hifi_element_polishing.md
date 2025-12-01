# Polishing Washu Pedigree with hifi and element data
## 1. PAN027

### 2. Comparing DV and DP:

https://t2t-consortium.slack.com/archives/C0696F9B6QK/p1730763498232509

ran as batch submission in: https://github.com/miramastoras/phoenix_batch_submissions/tree/main/workflows

```
Whole Genome:
unpolished:   54.141
Deepvariant:  54.2135
DeepPolisher: 54.1626

Confidence regions (~80% of genome):
unpolished:   70.981
Deepvariant:  71.5336
DeepPolisher: 71.2283
```

Results location:
```
/private/groups/patenlab/mira/phoenix_batch_executions/workflows/hprc_DeepPolisher
/private/groups/patenlab/mira/phoenix_batch_executions/workflows/hprc_polishing_QC
```
Testing hifi DV GQ filters: https://colab.research.google.com/drive/1LYvvQ8PtoPEb277LeXVxEMY1iglG52Ur?authuser=1#scrollTo=Ha9fwlb53hqy

Going with a GQ of 5
```
bcftools filter -Oz -i 'FORMAT/GQ[0]>5' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_patched/analysis/deepvariant_outputs/PAN027_patched_deepvariant.PASS.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_patched/analysis/deepvariant_outputs/PAN027_patched_deepvariant.GQ5.vcf.gz
```

### 3. Element Trio DeepVariant polishing with Q100 filters

Aligned all PAN027 reads to PAN027 maternal, and to paternal using bwa. Did the same for the parents (PAN010, PAN011)

https://github.com/miramastoras/phoenix_batch_submissions/tree/main/workflows/bwa

Ran deepvariant on each bam file

https://github.com/miramastoras/phoenix_batch_submissions/tree/main/workflows/deepvariant

replace sample name in gvcf
```
zcat /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_element_bwa_to_PAN027_mat/analysis/deepvariant_outputs/PAN027_element_bwa_to_PAN027_mat_deepvariant.g.vcf.gz | sed 's/default/PAN027/g' | bcftools view -Oz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_element_bwa_to_PAN027_mat/analysis/deepvariant_outputs/PAN027_element_bwa_to_PAN027_mat_deepvariant.fix.g.vcf.gz

zcat /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN027_mat/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN027_mat_deepvariant.g.vcf.gz | sed 's/default/PAN010/g' | bcftools view -Oz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN027_mat/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN027_mat_deepvariant.fix.g.vcf.gz

zcat /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN027_mat/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN027_mat_deepvariant.g.vcf.gz | sed 's/default/PAN011/g' | bcftools view -Oz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN027_mat/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN027_mat_deepvariant.fix.g.vcf.gz

zcat /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_element_bwa_to_PAN027_pat/analysis/deepvariant_outputs/PAN027_element_bwa_to_PAN027_pat_deepvariant.g.vcf.gz | sed 's/default/PAN027/g' | bcftools view -Oz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_element_bwa_to_PAN027_pat/analysis/deepvariant_outputs/PAN027_element_bwa_to_PAN027_pat_deepvariant.fix.g.vcf.gz

zcat /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN027_pat/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN027_pat_deepvariant.g.vcf.gz | sed 's/default/PAN010/g' | bcftools view -Oz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN027_pat/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN027_pat_deepvariant.fix.g.vcf.gz

zcat /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN027_pat/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN027_pat_deepvariant.g.vcf.gz | sed 's/default/PAN011/g' | bcftools view -Oz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN027_pat/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN027_pat_deepvariant.fix.g.vcf.gz
```

Run glnexus to merge the vcfs from each DV run
```
#!/bin/bash
#SBATCH --job-name=glnexus_PAN027
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=1:00:00

# maternal haplotype
time docker run \
  -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  ghcr.io/dnanexus-rnd/glnexus:v1.4.1 \
  glnexus_cli --config DeepVariant_unfiltered -t 32 \
  --dir /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/GLnexus.DB.mat \
  /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_element_bwa_to_PAN027_mat/analysis/deepvariant_outputs/PAN027_element_bwa_to_PAN027_mat_deepvariant.fix.g.vcf.gz \
  /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN027_mat/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN027_mat_deepvariant.fix.g.vcf.gz \
  /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN027_mat/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN027_mat_deepvariant.fix.g.vcf.gz \
  > /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_mat.element_deepvariant_1.6.1.glnexus.bcf

# paternal haplotype
time docker run \
  -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  ghcr.io/dnanexus-rnd/glnexus:v1.4.1 \
  glnexus_cli --config DeepVariant_unfiltered -t 32 \
  --dir /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/GLnexus.DB.pat \
  /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_element_bwa_to_PAN027_pat/analysis/deepvariant_outputs/PAN027_element_bwa_to_PAN027_pat_deepvariant.fix.g.vcf.gz \
  /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN027_pat/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN027_pat_deepvariant.fix.g.vcf.gz \
  /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN027_pat/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN027_pat_deepvariant.fix.g.vcf.gz \
  > /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_pat.element_deepvariant_1.6.1.glnexus.bcf
```
Convert to vcf
```
bcftools convert -Oz -o /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_pat.element_deepvariant_1.6.1.glnexus.vcf.gz /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_pat.element_deepvariant_1.6.1.glnexus.bcf

bcftools convert -Oz -o /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_mat.element_deepvariant_1.6.1.glnexus.vcf.gz /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_mat.element_deepvariant_1.6.1.glnexus.bcf
```
RTG mendelian to exclude any variants violating mendelian consistency, and also to phase the variant calls

```
rtg mendelian \
    -t /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/PAN027_raw_pat.fasta.sdf \
    -i /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_pat.element_deepvariant_1.6.1.glnexus.vcf.gz \
    --pedigree /private/groups/patenlab/mira/washu_pedigree_polishing/data/washu.ped \
    --output-consistent /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_pat.element_deepvariant_1.6.1.glnexus.phased.vcf.gz \
    --output-inconsistent /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_pat.element_deepvariant_1.6.1.glnexus.inconsistent.vcf.gz \
    --Xphase --all-records

rtg mendelian \
    -t /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/PAN027_raw_mat.fasta.sdf \
    -i /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_mat.element_deepvariant_1.6.1.glnexus.vcf.gz \
    --pedigree /private/groups/patenlab/mira/washu_pedigree_polishing/data/washu.ped \
    --output-consistent /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_mat.element_deepvariant_1.6.1.glnexus.phased.vcf.gz \
    --output-inconsistent /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_mat.element_deepvariant_1.6.1.glnexus.inconsistent.vcf.gz \
    --Xphase --all-records
```

Apply the Q100 filters:
- child GQ>20 and DP>20 and DP<double the mean coverage (90 in this case)
- phased to correct haplotype
- 1-2 bp indels

```
bcftools filter -i 'FORMAT/GQ[2]>20  && FORMAT/DP[2]>20 && MAX(FORMAT/DP[2])<180' /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_pat.element_deepvariant_1.6.1.glnexus.phased.vcf.gz | awk '$1 ~ /^#/ || ($12 ~ /^1/ || $12 ~ /^2/ )' | sed -r 's/1\|./1|1/; s/2\|./2|2/' | bcftools view -s PAN027 --trim-alt-alleles | bcftools norm -f /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/assembly.v1.0.PAN027.paternal.fa | awk '$1 ~ /^#/ || (length($4)+length($5)<5 && length($4)+length($5)>2)' > /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027_raw_pat.element_deepvariant_1.6.1.Q100_trio_filtered.final.vcf

bcftools filter -i 'FORMAT/GQ[2]>20  && FORMAT/DP[2]>20 && MAX(FORMAT/DP[2])<180' /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027-PAN010-PAN011.PAN027_raw_mat.element_deepvariant_1.6.1.glnexus.phased.vcf.gz | awk '$1 ~ /^#/ || ($12 ~ /^.\|1/ || $12 ~ /^.\|2/ )' | sed -r 's/.\|1/1|1/; s/.\|2/2|2/' | bcftools view -s PAN027 --trim-alt-alleles | bcftools norm -f /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/assembly.v1.0.PAN027.maternal.fa | awk '$1 ~ /^#/ || (length($4)+length($5)<5 && length($4)+length($5)>2)' > /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027_raw_mat.element_deepvariant_1.6.1.Q100_trio_filtered.final.vcf
```
Combine mat and pat vcfs
```
bcftools concat -Oz /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027_raw_mat.element_deepvariant_1.6.1.Q100_trio_filtered.final.vcf.gz /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027_raw_pat.element_deepvariant_1.6.1.Q100_trio_filtered.final.vcf.gz > /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027_raw_diploid.element_deepvariant_1.6.1.Q100_trio_filtered.final.vcf.gz
```
### Calculate overlaps between element and hifi final polishing vcfs

```
/private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027
```
Find intersections of Hifi and Element callsets. Require all alleles to be identical (-c none)
```
bcftools isec -c none -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_allele /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027_raw_diploid.element_deepvariant_1.6.1.Q100_trio_filtered.final.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_patched/analysis/deepvariant_outputs/PAN027_patched_deepvariant.GQ5.vcf.gz
```
Counts of intersections by allele
```
# 244 matching alleles
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_allele/0002.vcf | head -n 30
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_allele/0003.vcf | head -n 30

# 2558 just in element
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_allele/0000.vcf | head -n 30

# just in hifi 6231
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_allele/0001.vcf | head -n 30
```
Count of intersections by position
```
bcftools isec -c all -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_pos /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027_raw_diploid.element_deepvariant_1.6.1.Q100_trio_filtered.final.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_patched/analysis/deepvariant_outputs/PAN027_patched_deepvariant.GQ5.vcf.gz

# 248 matching alleles
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_pos/0002.vcf | head -n 30
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_pos/0003.vcf | head -n 30
```
Look at the 4 discordant calls in IGV
```
bedtools subtract -a /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_pos/0002.vcf -b  /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_allele/0002.vcf

bedtools subtract -a /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_pos/0003.vcf -b  /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_allele/0003.vcf
```
Subset bam files to 1000 bp around each site
```
bedtools subtract -a /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_pos/0003.vcf -b  /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/hifi_ele_isec_by_allele/0003.vcf | cut -f 1-2 | awk '{print $1, $2-10000, $2+10000}' OFS='\t' > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/discordant_alleles_slop_10kb.bed
```

Subset hifi PAN027 bam
```
samtools view \
  -L /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/discordant_alleles_slop_10kb.bed \
  -bh /private/groups/patenlab/mira/phoenix_batch_executions/workflows/hprc_DeepPolisher/PAN027_patched_verkko_model2/analysis/hprc_DeepPolisher_outputs/PAN027_patched_verkko_model2.hifi.to.diploid.asm.PHARAOH.bam \
  > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/PAN027_HiFi_PHARAOH.discordant.bam
```
Subset PAN027 element bam
```
samtools view \
  -L /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/discordant_alleles_slop_10kb.bed \
  -bh /private/groups/patenlab/mira/phoenix_batch_executions/workflows/bwa/PAN027_ele_to_PAN027_raw_mat/analysis/bwa_outputs/PAN027_ele_to_PAN027_raw_mat_bwa_0.7.17-r1188.sorted.bam \
  > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/PAN027_ele_to_PAN027_raw_mat_bwa.discordant.bam

samtools view \
  -L /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/discordant_alleles_slop_10kb.bed \
  -bh /private/groups/patenlab/mira/phoenix_batch_executions/workflows/bwa/PAN027_ele_to_PAN027_raw_pat/analysis/bwa_outputs/PAN027_ele_to_PAN027_raw_pat_bwa_0.7.17-r1188.sorted.bam \
  > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/PAN027_ele_to_PAN027_raw_pat_bwa.discordant.bam
```
Subset PAN010 element bams
```
samtools view \
  -L /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/discordant_alleles_slop_10kb.bed \
  -bh /private/groups/patenlab/mira/phoenix_batch_executions/workflows/bwa/PAN010_ele_to_PAN027_raw_mat/analysis/bwa_outputs/PAN010_ele_to_PAN027_raw_mat_bwa_0.7.17-r1188.sorted.bam \
  > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/PAN010_ele_to_PAN027_raw_mat_bwa.discordant.bam

s
samtools view \
  -L /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/discordant_alleles_slop_10kb.bed \
  -bh /private/groups/patenlab/mira/phoenix_batch_executions/workflows/bwa/PAN010_ele_to_PAN027_raw_pat/analysis/bwa_outputs/PAN010_ele_to_PAN027_raw_pat_bwa_0.7.17-r1188.sorted.bam \
  > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/PAN010_ele_to_PAN027_raw_pat_bwa.discordant.bam
```
Subset PAN011 bams
```
samtools view \
  -L /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/discordant_alleles_slop_10kb.bed \
  -bh /private/groups/patenlab/mira/phoenix_batch_executions/workflows/bwa/PAN011_ele_to_PAN027_raw_mat/analysis/bwa_outputs/PAN011_ele_to_PAN027_raw_mat_bwa_0.7.17-r1188.sorted.bam \
  > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/PAN011_ele_to_PAN027_raw_mat_bwa.discordant.bam

s
samtools view \
  -L /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/discordant_alleles_slop_10kb.bed \
  -bh /private/groups/patenlab/mira/phoenix_batch_executions/workflows/bwa/PAN011_ele_to_PAN027_raw_pat/analysis/bwa_outputs/PAN011_ele_to_PAN027_raw_pat_bwa_0.7.17-r1188.sorted.bam \
  > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/investigate_discordant_alleles/PAN011_ele_to_PAN027_raw_pat_bwa.discordant.bam
```
Slide deck: https://docs.google.com/presentation/d/1Gj3CFQmDMJQFKG8sJNNnPByPtQdSslFKEV7XD71Z9Og/edit#slide=id.g330dcff88f4_0_147

Removing discordant edits from the hifi callset, to move forward with the element calls
```
zcat /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_patched/analysis/deepvariant_outputs/PAN027_patched_deepvariant.GQ5.vcf.gz | grep "^#" > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_patched/analysis/deepvariant_outputs/PAN027_patched_deepvariant.GQ5.removed_discordant.vcf

zcat /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_patched/analysis/deepvariant_outputs/PAN027_patched_deepvariant.GQ5.vcf.gz | grep -v "86089335" | grep -v "113023980" | grep -v "44167493" | grep -v "3908681" | sed 's/default/PAN027/g' > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_patched/analysis/deepvariant_outputs/PAN027_patched_deepvariant.GQ5.removed_discordant.vcf
```
Combine hifi and element call set
```
bcftools concat -a -Oz /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027_raw_diploid.element_deepvariant_1.6.1.Q100_trio_filtered.final.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_patched/analysis/deepvariant_outputs/PAN027_patched_deepvariant.GQ5.removed_discordant.vcf.gz > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/PAN027_HiFi_DV_GQ5_element_DV_trio_Q100filt.final.vcf.gz

Polish
bcftools consensus -H2 -f /private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/diploid/assembly.v1.0.PAN027.diploid.fa /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/PAN027_HiFi_DV_GQ5_element_DV_trio_Q100filt.final.vcf.gz > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN027/polished.fasta
```

Intersect calls with Monika's homopolymer track
```
bedtools intersect -a /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN027_patched/analysis/deepvariant_outputs/PAN027_patched_deepvariant.GQ5.vcf.gz -b /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/assembly.v1.0.PAN027.diploid.mono+di.trf.minL10.slop5.bed | sort | uniq | wc -l

bedtools intersect -a /private/groups/patenlab/mira/washu_pedigree_polishing/glnexus/PAN027_raw_DeepVariant_bwa_element_unfiltered/PAN027_raw_diploid.element_deepvariant_1.6.1.Q100_trio_filtered.final.vcf.gz -b /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/assembly.v1.0.PAN027.diploid.mono+di.trf.minL10.slop5.bed | sort | uniq | wc -l
```

## 2. PAN028

Apply GQ filters to HiFi DV callset
```
bcftools view -Oz -f "PASS" /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN028_final_hifi_pharaoh_deepvariant.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN028_final_hifi_pharaoh_deepvariant.PASS.vcf.gz

bcftools filter -Oz -i 'FORMAT/GQ[0]>5' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN028_final_hifi_pharaoh_deepvariant.PASS.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN028_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz

# 11556 edits
```

Apply filters to element callset and extract only homozygous alt variants
```
# 9851 homalt edits
# restrict to 1-2 bp indels

bcftools filter -i 'FORMAT/GQ[0]>20  && FORMAT/DP[0]>20 && MAX(FORMAT/DP[0])<180' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_element_bwa_to_PAN028_hap1/analysis/deepvariant_outputs/PAN028_element_bwa_to_PAN028_hap1_deepvariant.vcf.gz | bcftools view -f "PASS" | awk '$1 ~ /^#/ || (length($4)+length($5)<5 && length($4)+length($5)>2)' | awk '/1\/1/ || /^#/ {print $0}' > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_element_bwa_to_PAN028_hap1/analysis/deepvariant_outputs/PAN028_element_bwa_to_PAN028_hap1_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf
# 9652

bcftools filter -i 'FORMAT/GQ[0]>20  && FORMAT/DP[0]>20 && MAX(FORMAT/DP[0])<180' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_element_bwa_to_PAN028_hap2/analysis/deepvariant_outputs/PAN028_element_bwa_to_PAN028_hap2_deepvariant.vcf.gz | bcftools view -f "PASS" | awk '$1 ~ /^#/ || (length($4)+length($5)<5 && length($4)+length($5)>2)' | awk '/1\/1/ || /^#/ {print $0}' > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_element_bwa_to_PAN028_hap2/analysis/deepvariant_outputs/PAN028_element_bwa_to_PAN028_hap2_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf
# 9723

# Combine element callsets
bcftools concat -Oz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_element_bwa_to_PAN028_hap2/analysis/deepvariant_outputs/PAN028_element_bwa_to_PAN028_hap2_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_element_bwa_to_PAN028_hap1/analysis/deepvariant_outputs/PAN028_element_bwa_to_PAN028_hap1_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf.gz > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/PAN028_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz
```
Intersect HiFi and element callsets
```
/private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028
```
```
bcftools isec -c none -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/hifi_ele_isec_by_allele /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/PAN028_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN028_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz
```
Counts of intersections by allele
```
# 974 matching alleles
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/hifi_ele_isec_by_allele/0002.vcf | head -n 30
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/hifi_ele_isec_by_allele/0003.vcf | head -n 30

# 18401 just in element
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/hifi_ele_isec_by_allele/0000.vcf | head -n 30

# just in hifi 10582
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/hifi_ele_isec_by_allele/0001.vcf | head -n 30
```
Count of intersections by position
```
bcftools isec -c all -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/hifi_ele_isec_by_pos /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/PAN028_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN028_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz

# 979 matching alleles
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/hifi_ele_isec_by_pos/0002.vcf | head -n 30
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/hifi_ele_isec_by_pos/0003.vcf | head -n 30
```
Combine HiFi and element callset. Remove any hifi edits that overlap element ones
```
# remove any hifi edits overlapping element ones
bedtools subtract -header \
    -A -a /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN028_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz \
    -b /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/PAN028_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN028_final_hifi_pharaoh_deepvariant.GQ5.ele_overlaps_removed.vcf

# 10563 edits

# merge vcf files
bcftools concat -a -Oz /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/PAN028_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN028_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN028_final_hifi_pharaoh_deepvariant.GQ5.ele_overlaps_removed.vcf.gz > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/PAN028_HiFi_DV_GQ5_element_DV_homalt_filtered.final.vcf.gz
# 29938 total edits
```

## 3. PAN010

Apply GQ filters to HiFi DV callset
```
bcftools view -Oz -f "PASS" /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN010_final_hifi_pharaoh_deepvariant.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN010_final_hifi_pharaoh_deepvariant.PASS.vcf.gz

bcftools filter -Oz -i 'FORMAT/GQ[0]>5' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN010_final_hifi_pharaoh_deepvariant.PASS.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN010_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz

# 9791 edits
```

Apply filters to element callset and extract only homozygous alt variants
```
bcftools filter -i 'FORMAT/GQ[0]>20  && FORMAT/DP[0]>20 && MAX(FORMAT/DP[0])<180' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN010_hap1/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN010_hap1_deepvariant.vcf.gz | bcftools view -f "PASS" | awk '$1 ~ /^#/ || (length($4)+length($5)<5 && length($4)+length($5)>2)' | awk '/1\/1/ || /^#/ {print $0}' > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN010_hap1/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN010_hap1_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf
# 5518

bcftools filter -i 'FORMAT/GQ[0]>20  && FORMAT/DP[0]>20 && MAX(FORMAT/DP[0])<180' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN010_hap2/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN010_hap2_deepvariant.vcf.gz | bcftools view -f "PASS" | awk '$1 ~ /^#/ || (length($4)+length($5)<5 && length($4)+length($5)>2)' | awk '/1\/1/ || /^#/ {print $0}' > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN010_hap2/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN010_hap2_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf
# 5551

mkdir -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/

# Combine element callsets
bcftools concat -Oz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN010_hap2/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN010_hap2_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_element_bwa_to_PAN010_hap1/analysis/deepvariant_outputs/PAN010_element_bwa_to_PAN010_hap1_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf.gz > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/PAN010_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz
```
Intersect HiFi and element callsets
```
bcftools isec -c none -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/hifi_ele_isec_by_allele /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/PAN010_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN010_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz
```
Counts of intersections by allele
```
# 743 matching alleles
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/hifi_ele_isec_by_allele/0002.vcf | head -n 30
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/hifi_ele_isec_by_allele/0003.vcf | head -n 30

# 10326 just in element
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/hifi_ele_isec_by_allele/0000.vcf | head -n 30

# just in hifi 9048
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/hifi_ele_isec_by_allele/0001.vcf | head -n 30
```
Count of intersections by position
```
bcftools isec -c all -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/hifi_ele_isec_by_pos /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/PAN010_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN010_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz

# 749 matching alleles
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/hifi_ele_isec_by_pos/0002.vcf | head -n 30
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/hifi_ele_isec_by_pos/0003.vcf | head -n 30
```
Combine HiFi and element callset. Remove any hifi edits that overlap element ones
```
# remove any hifi edits overlapping element ones
bedtools subtract -header \
    -A -a /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN010_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz \
    -b /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/PAN010_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN010_final_hifi_pharaoh_deepvariant.GQ5.ele_overlaps_removed.vcf

# merge vcf files
bcftools concat -a -Oz /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/PAN010_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN010_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN010_final_hifi_pharaoh_deepvariant.GQ5.ele_overlaps_removed.vcf.gz > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/PAN010_HiFi_DV_GQ5_element_DV_homalt_filtered.final.vcf.gz
# 20107 total edits
```
## 3. PAN011

Apply GQ filters to HiFi DV callset
```
bcftools view -Oz -f "PASS" /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.PASS.vcf.gz

bcftools filter -Oz -i 'FORMAT/GQ[0]>5' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.PASS.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz

# 13835 edits
```

Apply filters to element callset and extract only homozygous alt variants
```
bcftools filter -i 'FORMAT/GQ[0]>20  && FORMAT/DP[0]>20 && MAX(FORMAT/DP[0])<180' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN011_hap1/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN011_hap1_deepvariant.vcf.gz | bcftools view -f "PASS" | awk '$1 ~ /^#/ || (length($4)+length($5)<5 && length($4)+length($5)>2)' | awk '/1\/1/ || /^#/ {print $0}' > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN011_hap1/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN011_hap1_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf
# 18171

bcftools filter -i 'FORMAT/GQ[0]>20  && FORMAT/DP[0]>20 && MAX(FORMAT/DP[0])<180' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN011_hap2/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN011_hap2_deepvariant.vcf.gz | bcftools view -f "PASS" | awk '$1 ~ /^#/ || (length($4)+length($5)<5 && length($4)+length($5)>2)' | awk '/1\/1/ || /^#/ {print $0}' > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN011_hap2/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN011_hap2_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf
# 16222

mkdir -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/

# Combine element callsets
bcftools concat -Oz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN011_hap2/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN011_hap2_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_element_bwa_to_PAN011_hap1/analysis/deepvariant_outputs/PAN011_element_bwa_to_PAN011_hap1_deepvariant.PASS.GQ20.homalt.1_2bp_indels.vcf.gz > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz
```
Remove X and Y edits
```
# remove X and Y edits
zcat /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz | grep "^#" > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.noXY.vcf
zcat /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz | grep -v "^#" | grep -v "chrX" | grep -v "chrY" >> /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.noXY.vcf
```
Intersect HiFi and element callsets
```
bcftools isec -c none -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_allele /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.noXY.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz
```
Counts of intersections by allele
```
# 1228 matching alleles
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_allele/0002.vcf | head -n 30
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_allele/0003.vcf | head -n 30

# 31029 just in element
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_allele/0000.vcf | head -n 30

# just in hifi 12607
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_allele/0001.vcf | head -n 30
```
Count of intersections by position
```
bcftools isec -c all -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_pos /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.noXY.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz

# 1240 matching alleles
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_pos/0002.vcf | head -n 30
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_pos/0003.vcf | head -n 30
```
Combine HiFi and element callset. Remove any hifi edits that overlap element ones
```
# 46,696 total edit

# combine again
# remove any hifi edits overlapping element ones
bedtools subtract -header \
    -A -a /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz \
    -b /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.noXY.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.ele_overlaps_removed.vcf

# 12303 edits

# merge vcf files
bcftools concat -a -Oz /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.noXY.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.ele_overlaps_removed.vcf.gz > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_HiFi_DV_GQ5_element_DV_homalt_filtered.final.vcf.gz
# 44560 edits
```

## Look at distribution of element edits on each assembly

```
zcat /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq -c
zcat /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN010/PAN010_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq -c
zcat /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN028/PAN028_element_bwa_PASS.GQ20.homalt.1_2bp_indels.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq -c
```
Spot check in IGV
```
samtools view -bh \
  /private/groups/patenlab/mira/phoenix_batch_executions/workflows/hprc_DeepPolisher/PAN011_final/analysis/hprc_DeepPolisher_outputs/PAN011_final.hifi.to.diploid.asm.PHARAOH.bam \
  PAN011.chr19.haplotype1 > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/hprc_DeepPolisher/PAN011_final/analysis/hprc_DeepPolisher_outputs/PAN011_final.hifi.to.diploid.asm.PHARAOH.chr19_hap1.bam

# element
samtools view -bh \
  /private/groups/patenlab/mira/phoenix_batch_executions/workflows/bwa/PAN011_ele_to_PAN011_raw_hap1/analysis/bwa_outputs/PAN011_ele_to_PAN011_raw_hap1_bwa_0.7.17-r1188.sorted.bam \
  PAN011.chr19.haplotype1 > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/bwa/PAN011_ele_to_PAN011_raw_hap1/analysis/bwa_outputs/PAN011_ele_to_PAN011_raw_hap1_bwa_0.7.17-r1188.sorted.chr19.bam
```
https://docs.google.com/presentation/d/1QB-WbfH-v_EglPPjEoKEAgcRUtOWVacp-Mv2VJbhQXo/edit#slide=id.g3323fd9d6c7_0_32

### PAN011 X and Y polishing with element

Create new assembly with X and Y chromosome included with both haplotypes
```
cp /private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/diploid/assembly.v1.0.PAN011.haplotype2.fa .
cp /private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/diploid/assembly.v1.0.PAN011.haplotype1.fa .

cp /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/assembly.v1.0.PAN011.haplotype1.fa /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/assembly.v1.0.PAN011.haplotype1.XY.fa

samtools faidx /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/assembly.v1.0.PAN011.haplotype2.fa PAN011.chrY.haplotype2 >> /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/assembly.v1.0.PAN011.haplotype1.XY.fa

cp /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/assembly.v1.0.PAN011.haplotype2.fa /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/assembly.v1.0.PAN011.haplotype2.XY.fa

samtools faidx /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/assembly.v1.0.PAN011.haplotype1.fa PAN011.chrX.haplotype1 >> /private/groups/patenlab/mira/washu_pedigree_polishing/assemblies/assembly.v1.0.PAN011.haplotype2.XY.fa
```
Extract XY edits from vcf file, apply GQ and DP filters and take only homalt edits.

```
bcftools filter -i 'FORMAT/GQ[0]>20  && FORMAT/DP[0]>20 && MAX(FORMAT/DP[0])<180' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_ele_to_PAN011_raw_hap1_XY/analysis/deepvariant_outputs/PAN011_ele_to_PAN011_raw_hap1_XY_deepvariant.vcf.gz | bcftools view -f "PASS" | awk '$1 ~ /^#/ || (length($4)+length($5)<5 && length($4)+length($5)>2)' | awk '/1\/1/ || /^#/ {print $0}' | awk '/PAN011\.chrY\.haplotype2/ || '/PAN011\.chrX\.haplotype1/' || /^#/ {print $0}' > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_ele_to_PAN011_raw_hap1_XY/analysis/deepvariant_outputs/PAN011_ele_to_PAN011_raw_hap1_XY_deepvariant.vcf.PASS.GQ20.homalt.1_2bp_indels.vcf

# 260 Y edits
# 2013 X edits
```
Remove edits from the PAR regions

> coordinates share from Monika: https://docs.google.com/spreadsheets/d/1IsRRR3mOoVkKVdtNTjSAYz4jK5XZeD7VgK-uhchkK6U/edit?gid=2086031890#gid=2086031890

```
PAN011.chrX.haplotype1 0 2386588
PAN011.chrX.haplotype1 153719008 154053196
PAN011.chrY.haplotype2 0 2427643
PAN011.chrY.haplotype2 46873975 47210280
```
```sh
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups pegi3s/bedtools bedtools subtract -header -a /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_ele_to_PAN011_raw_hap1_XY/analysis/deepvariant_outputs/PAN011_ele_to_PAN011_raw_hap1_XY_deepvariant.vcf.PASS.GQ20.homalt.1_2bp_indels.vcf -b /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_ele_to_PAN011_raw_hap1_XY/analysis/deepvariant_outputs/PAR.bed > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_ele_to_PAN011_raw_hap1_XY/analysis/deepvariant_outputs/PAN011_ele_to_PAN011_raw_hap1_XY_deepvariant.vcf.PASS.GQ20.homalt.1_2bp_indels.exclPAR.vcf

# excludes 127 edits

# 1956 X edits
# 190 Y edits
```

Recombine with final element callset
```sh
bcftools concat -Oz -a /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_ele_to_PAN011_raw_hap1_XY/analysis/deepvariant_outputs/PAN011_ele_to_PAN011_raw_hap1_XY_deepvariant.vcf.PASS.GQ20.homalt.1_2bp_indels.exclPAR.vcf.gz /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.noXY.vcf.gz > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.with_XY.vcf.gz
# 34403
```

Re-run merging with hifi edits
```sh
# remove any hifi edits overlapping element ones
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups pegi3s/bedtools bedtools subtract -header \
    -A -a /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz \
    -b /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.with_XY.vcf.gz > /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.ele_overlaps_removed.with_ele_XY.vcf

# 12322 edits

# merge vcf files
bcftools concat \
    -a -Oz /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.with_XY.vcf.gz \
    /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.ele_overlaps_removed.with_ele_XY.vcf.gz \
    > /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_HiFi_DV_GQ5_element_DV_homalt_filtered.final.withXY.vcf.gz

# 46725 edits
```

Intersect HiFi and element callsets
```
bcftools isec -c none -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_allele /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.with_XY.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz
```
Counts of intersections by allele
```
# 1489 matching alleles
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_allele/0002.vcf | head -n 30
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_allele/0003.vcf | head -n 30

# 32914 just in element
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_allele/0000.vcf | head -n 30

# just in hifi 12346
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_allele/0001.vcf | head -n 30
```
Count of intersections by position
```
bcftools isec -c all -p /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_pos /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/PAN011_element_bwa_PASS.GQ20.homalt.1_2bp_indels.with_XY.vcf.gz /private/groups/patenlab/mira/phoenix_batch_executions/workflows/deepvariant/PAN011_final_hifi_pharaoh/analysis/deepvariant_outputs/PAN011_final_hifi_pharaoh_deepvariant.GQ5.vcf.gz

# 1505 matching alleles
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_pos/0002.vcf | head -n 30
bcftools stats /private/groups/patenlab/mira/washu_pedigree_polishing/analysis/intersect_hifi_element_PAN011/hifi_ele_isec_by_pos/0003.vcf | head -n 30

1505-1489
```
### Adding additional kmer QV calculations to the manuscript

Concatenating the merqury QV results
```sh
cd /private/groups/patenlab/mira/washu_pedigree_polishing/QV_compare_technologies

cat /private/groups/patenlab/mira/washu_pedigree_polishing/QV_compare_technologies/QV_tech_list.txt | while read line ; do
  hap1=`sed -n '1p' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/merqury/${line}/analysis/merqury_outputs/*.polished.merqury.qv | cut -f4`
  hap2=`sed -n '2p' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/merqury/${line}/analysis/merqury_outputs/*.polished.merqury.qv | cut -f4`
  both=`sed -n '3p' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/merqury/${line}/analysis/merqury_outputs/*.polished.merqury.qv | cut -f4`
  echo ${line},${hap1},${hap2},${both}>> /private/groups/patenlab/mira/washu_pedigree_polishing/QV_compare_technologies/polished_merqury.csv
done

# unpolished
cat /private/groups/patenlab/mira/washu_pedigree_polishing/QV_compare_technologies/QV_tech_list_raw.txt | while read line ; do
  hap1=`sed -n '1p' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/merqury/${line}/analysis/merqury_outputs/*.merqury.qv | cut -f4`
  hap2=`sed -n '2p' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/merqury/${line}/analysis/merqury_outputs/*.merqury.qv | cut -f4`
  both=`sed -n '3p' /private/groups/patenlab/mira/phoenix_batch_executions/workflows/merqury/${line}/analysis/merqury_outputs/*.merqury.qv | cut -f4`
  echo ${line},${hap1},${hap2},${both}>> /private/groups/patenlab/mira/washu_pedigree_polishing/QV_compare_technologies/raw_merqury.csv
done
```
