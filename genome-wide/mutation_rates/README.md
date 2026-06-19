# Mutation Rates Analysis — PBS HPC Port

**Notebooks:** `step2c_mutation_rates.ipynb`, `validate_mutation_rates.ipynb`, `diagnose_length_mismatch.ipynb`, `2b_merge_variants.ipynb`
**CLI:** `mutation_rates_cli.py`, `preprocess_for_mutation_rates.py`

This codebase estimates mutation rates in the **PAN027** assembly by stratifying variants across high-confidence (easy) and difficult (hard) regions, as well as assembly-specific lower-confidence regions.

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

---

## Directory Layout

| File | Purpose |
|------|---------|
| `2b_merge_variants.ipynb` | Merges per-block VCFs from switch error pipeline's step 2 into `merged_variants.vcf.gz` and creates `all_blocks.bed`. |
| `2c_mutation_rates.ipynb` | Main mutation-rates calculation notebook. Stratifies variants by haplotype, region scope, and variant type. |

---

## Required Inputs

All inputs are expected to be available in the working directory unless otherwise specified.

### 1. Variant data
| File | Description |
|------|-------------|
| `merged_variants.vcf.gz` + `.tbi` | Merged, bgzipped, tabix-indexed variant calls from step 2. Should exclude sites from problematic grandparent or maternal regions. |

### 2. Transmitted-block coverage
| File | Description |
|------|-------------|
| `all_blocks.bed` | BED file describing all genomic regions covered by transmitted blocks in the pedigree (derived from `results.json` created in the switch_errors pipeline). |

### 3. Switch-error regions annotation
| File | Description |
|------|-------------|
| `sw_prone_regions.bed` | BED-formatted regions flagged as switch-error–prone by the `6_summarize.py` step of switch_errors pipeline. |

### 4. PAN027 problematic regions
| File | Description |
|------|-------------|
| `flagger.PAN027.ONT.bed` | Low-confidence regions identified by Flagger (ONT data). |
| `flagger.PAN027.hifi.bed` | Low-confidence regions identified by Flagger (HiFi data). |
| `problematic.PAN027.bed` | Combined Low-confidence regions specific to PAN027 assembly. |

### 5. Unreliable grandparental variants (optional)
| File | Description |
|------|-------------|
| `problematic_grandparents.vcf.gz` + `.tbi` | Variants originating from low-confidence regions in the grandparents' assemblies.|

### 6. GIAB-derived stratifications (lifted to PAN027)
These are created by mapping GIAB CHM13 stratifications onto PAN027 haplotypes.

**Maternal haplotype:**
| File | Description |
|------|-------------|
| `chm13_to_PAN027_mat.easy.PAN027names.bed` | Easy/high-confidence regions |
| `chm13_to_PAN027_mat.hard.PAN027names.bed` | Difficult regions (complement of easy) |

**Paternal haplotype:**
| File | Description |
|------|-------------|
| `chm13_to_PAN027_pat.easy.PAN027names.bed` | Easy/high-confidence regions |
| `chm13_to_PAN027_pat.hard.PAN027names.bed` | Difficult regions (complement of easy) |

**Global (union of maternal + paternal):**
| File | Description |
|------|-------------|
| `PAN027.v1.1.easy.bed` | Union of maternal + paternal easy regions |
| `PAN027.v1.1.difficult.bed` | Union of maternal + paternal difficult regions |

**GIAB source reference:** 
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/CHM13@all/Union/CHM13_notinalldifficultregions.bed.gz

### 7. Padding parameter
```python
PAD = 0  # Amount of padding added to flagged/problematic PAN027 regions
```

---

## Outputs

| File | Description |
|------|-------------|
| `mutation_rates_summary.tsv` | Coarse (whole-genome) mutation rates stratified by haplotype × region × variant type. |
| `mutation_rates_per_chromosome.tsv` | Fine (per-chromosome) mutation rates for all chromosomes. |
| `mutation_rates_chromosomes_all.tsv` | Reshaped per-chromosome table (all regions). |
| `mutation_rates_chromosomes_easy.tsv` | Reshaped per-chromosome table (easy regions only). |
| `mutation_rates_chromosomes_hard.tsv` | Reshaped per-chromosome table (hard regions only). |
| `validation_*.tsv` | Validation results comparing output TSVs against source BEDs. |

---

## Column Definitions

### Length columns
| Column | Description |
|--------|-------------|
| `transmitted_length` | Total length of transmitted blocks in the given scope (unfiltered denominator). Does not exclude problematic regions. |
| `polished_length` | Callable/non-problematic region size after filtering out switch-error–prone regions, Flagger regions, and padded problematic regions. **Used as the mutation-rate denominator.** |
| `region_scope_len` | Alias for `polished_length` in some outputs; represents the effective callable length for the given haplotype × region combination. |

### Variant count columns
| Column | Description |
|--------|-------------|
| `variants_count` | Raw count of variants in the given scope (before subtracting unreliable grandparent variants). |
| `variants_count_subtracted` | Count after subtracting unreliable grandparent variants. |
| `variants_per_Mb` | Mutation rate: `variants_count_subtracted / (polished_length / 1_000_000)`. |

### Stratification columns
| Column | Values | Description |
|--------|--------|-------------|
| `haplotype_scope` | `all`, `maternal`, `paternal` | Which haplotype's transmitted blocks are considered. |
| `region_scope` | `all`, `easy`, `hard` | GIAB-based region stratification. |
| `variant_type` | `all`, `SNVs`, `1bp_indels`, `>1bp_indels` | Variant type classification. Indels split by length difference: `abs(strlen(ALT)-strlen(REF))`. |

