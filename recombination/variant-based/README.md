## Prerequisites for the variant-calling based approach:

## 1. Align HiFi and ONT reads to haplotype assemblies

Reads are aligned separately to the **paternal** and **maternal** assemblies using `minimap2`.

### HiFi → Paternal assembly

```
minimap2 -ax map-hifi -t24 assembly.v1.0.PAN027.paternal.fa \
  HG06804_MGISTL-PAN011_grandfather/PacBio_HiFi/m64136_210506_183715.hifi_reads.fastq.gz \
  HG06804_MGISTL-PAN011_grandfather/PacBio_HiFi/m64136_210508_020125.hifi_reads.fastq.gz \
  HG06804_MGISTL-PAN011_grandfather/PacBio_HiFi/m64136_210509_092431.hifi_reads.fastq.gz \
  HG06804_MGISTL-PAN011_grandfather/PacBio_HiFi/m64136_210510_164711.hifi_reads.fastq.gz \
  -o PAN011-hifi_to_PAN027-assembly.v1.0.paternal.sam
```
### HiFi → Maternal assembly

```
minimap2 -ax map-hifi -t24 assembly.v1.0.PAN027.maternal.fa \
  HG06803_MGISTL-PAN010_grandmother/PacBio_HiFi/m64043_210416_171718.hifi_reads.fastq.gz \
  HG06803_MGISTL-PAN010_grandmother/PacBio_HiFi/m64043_210420_151018.hifi_reads.fastq.gz \
  HG06803_MGISTL-PAN010_grandmother/PacBio_HiFi/m64043_210418_021436.hifi_reads.fastq.gz \
  HG06803_MGISTL-PAN010_grandmother/PacBio_HiFi/m64043_210421_224017.hifi_reads.fastq.gz \
  -o PAN010-hifi_to_PAN027-assembly.v1.0.maternal.sam
```
### ONT → Paternal assembly

```
minimap2 -ax map-ont -t24 assembly.v1.0.PAN027.paternal.fa \
  HG06804_MGISTL-PAN011_grandfather/nanopore/03_15_22_R941_HG06804_1_Guppy_5.0.13_sup.fastq.gz \
  HG06804_MGISTL-PAN011_grandfather/nanopore/03_15_22_R941_HG06804_2_Guppy_5.0.13_sup.fastq.gz \
  HG06804_MGISTL-PAN011_grandfather/nanopore/03_15_22_R941_HG06804_3_Guppy_5.0.13_sup.fastq.gz \
  HG06804_MGISTL-PAN011_grandfather/nanopore/03_15_22_R941_HG06804_4_Guppy_5.0.13_sup.fastq.gz \
  HG06804_MGISTL-PAN011_grandfather/nanopore/03_15_22_R941_HG06804_5_Guppy_5.0.13_sup.fastq.gz \
  HG06804_MGISTL-PAN011_grandfather/nanopore/03_15_22_R941_HG06804_6_Guppy_5.0.13_sup.fastq.gz \
  HG06804_MGISTL-PAN011_grandfather/nanopore/03_15_22_R941_HG06804_7_Guppy_5.0.13_sup.fastq.gz \
  -o PAN011-ONT_to_PAN027-assembly.v1.0.paternal.sam
```
### ONT → Maternal assembly

```
minimap2 -ax map-ont -t24 assembly.v1.0.PAN027.maternal.fa \
  HG06803_MGISTL-PAN010_grandmother/nanopore/03_15_22_R941_HG06803_1_Guppy_5.0.13_sup.fastq.gz \
  HG06803_MGISTL-PAN010_grandmother/nanopore/03_15_22_R941_HG06803_4_Guppy_5.0.13_sup.fastq.gz \
  HG06803_MGISTL-PAN010_grandmother/nanopore/03_15_22_R941_HG06803_2_Guppy_5.0.13_sup.fastq.gz \
  HG06803_MGISTL-PAN010_grandmother/nanopore/03_15_22_R941_HG06803_5_Guppy_5.0.13_sup.fastq.gz \
  HG06803_MGISTL-PAN010_grandmother/nanopore/03_15_22_R941_HG06803_3_Guppy_5.0.13_sup.fastq.gz \
  HG06803_MGISTL-PAN010_grandmother/nanopore/03_15_22_R941_HG06803_6_Guppy_5.0.13_sup.fastq.gz \
  -o PAN010-ONT_to_PAN027-assembly.v1.0.maternal.sam
```
## 2. Call variants in the assemblies
```
singularity exec deepvariant_1.5.0.sif \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref assembly.v1.0.PAN027.paternal.fa \
  --num_shards 24 \
  --reads PAN011-hifi_to_PAN027-assembly.v1.0.paternal.bam \
  --sample_name PAN011-to-PAN027-assembly.v1.0-paternal \
  --output_vcf PAN011_to_PAN027-assembly.v.0-paternal_deepvariant.vcf \
  --intermediate_results_dir /temp-PAN011-to-PAN027-assembly.v1.0-paternal/

singularity exec deepvariant_1.5.0.sif \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref assembly.v1.0.PAN027.maternal.fa \
  --num_shards 24 \
  --reads PAN010-hifi_to_PAN027-assembly.v1.0.maternal.bam \
  --sample_name PAN010-to-PAN027-assembly.v1.0-maternal \
  --output_vcf PAN010_to_PAN027-assembly.v.0-maternal_deepvariant.vcf \
  --intermediate_results_dir /temp-PAN010-to-PAN027-assembly.v1.0-maternal/
```
## 3. Phase variants
```
whatshap phase --ignore-read-groups \
  -o PAN011_to_PAN027-assembly.v1.0-paternal_deepvariant_phased_ont-hifi.vcf \
  --reference=assembly.v1.0.PAN027.paternal.fa \
  PAN011_to_PAN027-assembly.v1.0-paternal_deepvariant.vcf \
  PAN011-ONT_to_PAN027-assembly.v1.0.paternal.bam \
  PAN011-hifi_to_PAN027-assembly.v1.0.paternal.bam

whatshap phase --ignore-read-groups \
  -o PAN010_to_PAN027-assembly.v1.0-maternal_deepvariant_phased_ont-hifi.vcf \
  --reference=assembly.v1.0.PAN027.maternal.fa \
  PAN010_to_PAN027-assembly.v1.0-maternal_deepvariant.vcf \
  PAN010-ONT_to_PAN027-assembly.v1.0.maternal.bam \
  PAN011-hifi_to_PAN027-assembly.v1.0.maternal.bam
```

## 4. Filter phased variants
#### Keep only phased variant calls
`bcftools view --phased --threads 2 \
  PAN011_to_PAN027-assembly.v1.0-paternal_deepvariant_phased_ont-hifi.vcf \
  -o PAN011_to_PAN027-assembly.v1.0-paternal_deepvariant_phased_ont-hifi_phasedcalls.vcf
`
#### Keep only calls with minimum QV of 30
`vcftools --gzvcf PAN011_to_PAN027-assembly.v1.0-paternal_deepvariant_phased_ont-hifi_phasedcalls.vcf \
  --minQ 30 \
  --out PAN011_to_PAN027-assembly.v1.0-paternal_deepvariant_phased_ont-hifi_phasedcalls_q30.vcf \
  --recode
`
#### Keep only calls with minimum depth value of 65
`vcftools --gzvcf PAN011_to_PAN027-assembly.v1.0-paternal_deepvariant_phased_ont-hifi_phasedcalls_q30.vcf.recode.vcf \
  --max-meanDP 65 \
  --out PAN011_to_PAN027-assembly.v1.0-paternal_deepvariant_phased_ont-hifi_phasedcalls_q30_maxdp65_site \
  --recode
`
## 5. Detect recombination breakpoints in VCF files:

`python find_recombination_vcf.py PAN010_to_PAN027-assembly.v1.0-maternal_deepvariant_phased_ont-hifi_phasedcalls_q30_maxdp65_site.recode.vcf recomb_PAN027.v1.0-maternal_q30_maxdp65_100kb.tsv > recombination/boundaries_input_PAN027.v1.0-maternal.tsv`

`python create_inputfile_recombination_PAN027.v1.0.py boundaries_input_PAN027.v1.0-maternal.tsv bedfiles_PAN027.v1.0-maternal/recomb_boundaries bedfiles_PAN027.v1.0-maternal/phaseset_boundaries`

`python add_offset_to_boundaries_PAN027.v1.0.py ragtag-PAN027-maternal/ragtag.scaffold.agp bedfiles_PAN027.v1.0-maternal`


