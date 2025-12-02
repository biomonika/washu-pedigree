# De-novo Mutation Enrichment Analysis  
**Notebook:** `enrichment.ipynb`

This Jupyter notebook evaluates whether variants produced by the **switch_errors** pipeline are enriched or depleted with respect to various genomic annotations. It measures the proximity of observed variants to annotation features, compares this against permutation-based expectations, and summarizes enrichment statistics and significance.

---

## Environment Setup

The notebook is designed to run inside a dedicated Conda environment.

### 1. Create the environment
```bash
conda env create -f environment.yml
```

### 2. Install Python package requirements (inside the Conda environment)
```bash
pip install -r requirements.txt
```

## Required Input Data

The notebook expects the following directory structure:

```
data/
  blocks/
    results.json
  variants/
    ... VCFs from step 2_call_variants.py ...
  experiment/
    annots_folder/
      *.bed   # annotation files for PAN027
```

### 1. Variants

Variants must come from step 2_call_variants.py of the switch-errors pipeline.
They must be available under:

`./data/variants/`

### 2. Switch error detection pipeline metadata

The results.json of the switch-errors pipeline must be located at:

`data/blocks/results.json`

### 3. Annotation BED files

Annotation files for PAN027 must be placed in:

`./data/experiment/annots_folder/`

Each BED file must contain a column used to define annotation types, determined by:

    ANNOTS_FIELD — the column index or key in the BED record used to extract the annotation value
    ANNOTS_TYPE_LAMBDA — a callable transforming that value into the normalized annotation class used in the notebook

Examples are provided within the notebook

## Analysis Workflow
### 1. Mean distance to annotations

For every annotation type, the notebook computes the distance from 
each observed variant to the nearest corresponding annotation feature

### 2. Permutation testing

Variants are permuted within their transmitted blocks and for each permutation, 
the mean distance to annotation features is recomputed.
The result is a distribution of expected mean distances per annotation type

### 3. Enrichment & statistical significance

For each annotation type, the notebook computes:

    Fold enrichment / depletion
    Empirical p-values
    Ranking of annotation categories by effect size and significance
    Plots summarizing the enrichment patterns

## Output Files

The notebook produces the following outputs in the working directory:

| File | Description                                                                                        |
|------|----------------------------------------------------------------------------------------------------|
| enrichment_fold_enrichment_bar.png     | Bar plot showing enrichment/depletion and statistical significance for each annotation type        |
| enrichment_topk_perm_histograms.png     | Histograms of permuted mean distances for the top enriched/depleted annotation types               |
| enrichment_detailed.tsv     | Per-annotation statistics: observed mean distance, expected distribution, p-value, fold enrichment |
| enrichment_perm_means.npy     | Raw permutation mean-distance data (numpy array)                                                   |
