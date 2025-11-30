# Introduction & Methods
This report focuses on detecting novel repeat insertions within syntenic regions shared across a three-generation pedigree. We analysed shared syntenic blocks between grandparents, mother, and granddaughter assemblies. Next, we called variants in shared blocks between generations. We filtered variants to only insertion longer than 50bp and determined whether these insertions represent novel TE copies by intersecting with repeat annotation created by RepeatMasker. This resulted in insertion candidates that were manually checked in IGV.

In mother -> granddaughter variant file, mother assembly was used as reference, and therefore could not be used to cross-check with repeat track from the granddaughter assembly. For this reason, variants file for the mother -> granddaughter comparison was first converted to use granddaughter assembly as ref using convert_vcf_query_to_ref.py.

# Code
generation 1 - grandparents to mother
```bash
file="diploid_generation1_shifted.20250901.vcf.gz.20250902"
# select long insertions
cat data/polished/${file}.vcf | awk 'length($4) > 50' > data/polished/${file}.long.vcf 
# add header to use with bedtools and view in igv
cat data/header.vcf data/polished/${file}.long.vcf > data/polished/${file}.long.header.vcf
# instersect long insertions with repeat annotations
bedtools intersect -a data/polished/PAN027_diploid.RepeatMasker.dedup.bed -b data/polished/${file}.long.header.vcf | grep -v "Satellite\|Simple_repeat\|Low_complexity\|ACRO" > data/polished/te_onsertion_candidates_generation1.bed
```

generation 2 - mother to granddaughter
```bash
cat data/polished/PAN028_*RepeatMasker.bed | sort -k1,1 -k2,2n > data/polished/PAN028_diploid.RepeatMasker.bed
python workflow/convert_vcf_query_to_ref.py   data/polished/diploid_generation2_shifted.20250901.vcf.gz.20250902.vcf   data/polished/diploid_generation2_PAN028.CLI.vcf   -b data/polished/maternal.blocks.formatted.filtered.tsv   -b data/polished/paternal.blocks.formatted.filtered.tsv
file="diploid_generation2_PAN028.CLI"
# select long insertions
cat data/polished/${file}.vcf | awk 'length($4) > 50' > data/polished/${file}.long.vcf 
# add header to use with bedtools and view in igv
cat data/header.vcf data/polished/${file}.long.vcf > data/polished/${file}.long.header.vcf
# instersect long insertions with repeat annotations
bedtools intersect -a data/polished/PAN028_diploid.RepeatMasker.bed -b data/polished/${file}.long.header.vcf | grep -v "Satellite\|Simple_repeat\|Low_complexity\|ACRO" > data/polished/te_onsertion_candidates_generation2.bed
```
# Results

By intersecting insertions in shared blocks with repeat track, we found 6 insertion candidates in generation 1 (grandparents -> mother) and 5 candidates in generation 2 (mother -> granddaughter). After manually checking each candidate in Integrative Genomics Viewer (IGV), all candidates were found to be a result of incorrect haplotype phasing (element was present in the second haplotype of the parent, which can be a result of an switch error or gene conversion event) or incorrect variant call.