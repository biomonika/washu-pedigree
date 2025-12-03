# Mutation Rates Analysis  
**Notebook:** `mutation_rates.ipynb`

This Jupyter notebook estimates mutation rates in the **PAN027** assembly by stratifying 
variants across high-confidence and difficult regions (based on GIAB definitions), 
as well as assembly-specific problematic regions. 

---

## Environment Setup

We provide conda environment for reproducibility:

### 1. Create Conda environment
```bash
conda env create -f environment.yml
```
### 2. Install Python dependencies
```bash
pip install -r requirements.txt
```

## Required Inputs

All inputs are expected to be available in the notebook’s working directory unless otherwise specified.

### 1. Variant data

    variants.vcf.gz
Concatenated, filtered variant calls from step 2_call_variants.py of the switch_errors pipeline.
These variants should exclude sites originating from problematic grandparent or maternal regions.

### 2. Switch-error regions annotation

    sw_prone_regions.bed
BED-formatted regions flagged as switch-error–prone by the 6_summarize.py step of switch_errors pipeline.

### 3. PAN027 problematic regions

    flagger.PAN027.ONT.bed
    flagger.PAN027.hifi.bed
    problematic.PAN027.bed

Problematic regions specific to the PAN027 assembly, identified by flagger and nucflag

### 4. Unreliable grandparental variants

    problematic_grandparents.vcf.gz
Variants originating from problematic regions in the grandparents’ assemblies.

### 5. GIAB-derived stratifications (lifted to PAN027)

These are created by mapping GIAB CHM13 stratifications onto PAN027 haplotypes.

Maternal haplotype:
    chm13_to_PAN027_mat.easy.mapped.PAN027names.bed – easy/high-confidence regions
    chm13_to_PAN027_mat.hard.PAN027names.bed – difficult regions (complement of easy)

Paternal haplotype:
    chm13_to_PAN027_pat.easy.mapped.PAN027names.bed – easy/high-confidence regions
    chm13_to_PAN027_pat.hard.PAN027names.bed – difficult regions (complement of easy)

GIAB source reference:
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/CHM13@all/Union/CHM13_notinalldifficultregions.bed.gz

### 6. Transmitted-block coverage

    all_blocks.bed
BED file describing all genomic regions covered by transmitted blocks in the pedigree (derived from results.json created in the switch_errors pipeline).

### 7. Padding parameter
    PAD = 20000

Amount of padding added to flagged/problematic PAN027 regions

## Analysis Overview

The notebook performs the following steps:
### 1. Load variant calls and mask unreliable regions

Exclude variants originating from:
        switch-error–prone regions (sw_prone_regions.bed)
        problematic PAN027 regions (ONT/HiFi/combined)
        unreliable grandparent regions

Apply padding (PAD) to problematic regions

### 2. Intersect variants with GIAB stratifications

Mutation counts are computed separately for:

    maternal easy
    maternal difficult
    paternal easy
    paternal difficult

Stratified mutation rates are computed by dividing counts by the callable/non-problematic region sizes in each category

### 3. Adjust for transmitted-block coverage

Only genomic regions covered by transmitted blocks (all_blocks.bed) are considered

### 4. Summarize mutation rates

The notebook reports mutation rates per (maternal regions, paternal regions, all regions) x (easy regions, difficult regions, all regions) x (SNVs, indels, all variants)