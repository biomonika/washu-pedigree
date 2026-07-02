# PAN027 Pangenome – De novo Mutation Analysis

## Samples

| Sample | Role |
|--------|------|
| PAN027 | Mother (index individual) |
| PAN010 | Grandmother (maternal line) |
| PAN011 | Grandfather (paternal line) |
| PAN028 | Granddaughter |
| CHM13  | T2T reference genome |

---

## Scripts

| Script | Description |
|--------|-------------|
| `run_cactus.sh` | Build pangenomes for all chromosomes (PAN027_mat + PAN027_pat projections, CHM13 included as sample) |
| `analyse_block.maternal.sh` | Filter de novo candidates – maternal line, template script called per chromosome |
| `analyse_block.paternal.sh` | Filter de novo candidates – paternal line, template script called per chromosome |
| `run_analyse_maternal.sh` | Run maternal filtering for all chromosomes |
| `run_analyse_paternal.sh` | Run paternal filtering for all chromosomes |
| `filter_vcf_into_easy_and_difficult.sh` | Split filtered VCFs into easy and difficult genomic regions |
| `make_stats.per_chr.py` | Compute mutation statistics per chromosome |
| `create_chrom_map.py` | Create chromosome map visualization (per chromosome) |
| `run_make_stats_all.sh` | Run statistics and maps for all chromosomes |
| `vg_pos_lookup.py` | Helper: maps VCF positions to alt haplotype coordinates via vg graph |

### Configuration

Each script contains a clearly marked `USER CONFIGURATION` block at the top. Edit these variables before running:

**`run_cactus.sh`**
```bash
PROJECT_DIR="/path/to/your/project"   # root project directory
SCRATCH_DIR="/path/to/your/scratch"   # scratch space (~100 GB needed)
IMAGE="quay.io/..."                   # Cactus Podman image version
DATA_DIR="${PROJECT_DIR}/data"        # FASTA input files root
```

**`analyse_block.maternal/paternal.sh`**
```bash
BLOCKS_DIR="/path/to/maternal_transmitted_blocks_per_chr"   # folder with per-chr block TSVs
CACTUS_RESULTS_DIR="/path/to/cactus/results"                # folder with PAN027.chrN.project_PAN027_mat/ subfolders
OUTPUT_DIR="/path/to/output"                                 # where filtered VCFs will be written
ANNOTATIONS_DIR="/path/to/annotations"                       # folder with BED annotation files
```

**`run_analyse_maternal/paternal.sh`**
```bash
PROJECT_DIR="/path/to/project"
OUT_BASE="${PROJECT_DIR}/results_final_filter_mother_gp"
STAGING_DIR="${PROJECT_DIR}/results_filtered.based_on_blocks.PAN027only"
```

**`filter_vcf_into_easy_and_difficult.sh`**
```bash
CONDA_SH="/path/to/conda.sh"                                 # path to conda.sh
CONDA_ENV="/path/to/conda/env"                               # conda environment path or name
INPUT_DIR="results_final_filter_mother_gp"                   # folder with filtered VCFs (input)
OUTPUT_DIR="results_final_filter_mother_gp_easy_difficult"   # output folder
EASY_BED="PAN027.v1.1.easy.bed"                             # easy regions BED file
DIFFICULT_BED="PAN027.v1.1.difficult.bed"                    # difficult regions BED file
```

**`make_stats.per_chr.py`**
```python
BLOCKS_DIR_MAT  = "/path/to/maternal_transmitted_blocks_per_chr"
BLOCKS_DIR_PAT  = "/path/to/paternal_transmitted_blocks_per_chr"
ANNOTATIONS_DIR = "/path/to/annotations"
```

**`create_chrom_map.py`**
```python
BLOCKS_DIR_MAT = "/path/to/maternal_transmitted_blocks_per_chr"
BLOCKS_DIR_PAT = "/path/to/paternal_transmitted_blocks_per_chr"
```

**`run_make_stats_all.sh`**
```bash
SCRIPTS_DIR="/path/to/scripts"          # folder where the Python scripts live
OUT_BASE="/path/to/results_folder"      # results folder to process
```

---

## Workflow

### Step 1 – Pangenome construction (`run_cactus.sh`)

Pangenomes were built per chromosome using **Cactus-pangenome v3.1.4** via Podman container.  
The script iterates over chromosomes 1–22 + X sequentially (one finishes before the next starts).  
CHM13 is included as a sample in the pangenome. For each chromosome, two pangenomes are built:
- projected onto **PAN027_mat** (maternal haplotype as reference)
- projected onto **PAN027_pat** (paternal haplotype as reference)

Each run produces: HAL alignment, GBZ graph, GFA, VCF (raw + wave-decomposed), PAF, GAF.  
Output: `results/PAN027.chr{N}.project_PAN027_{mat|pat}/`

Variants were called with `--vcf --vcfwave`. Wave decomposition (`vcfwave`) decomposes complex alleles into primitive SNVs and indels, writing `TYPE` and `LEN` fields into the VCF INFO column.

> **Note:** Set `PROJECT_DIR`, `SCRATCH_DIR`, and `DATA_DIR` in the USER CONFIGURATION block at the top of the script before running.

---

### Step 2 – De novo variant filtering (`run_analyse_maternal/paternal.sh`)

Run the appropriate script to process all chromosomes:

```bash
./run_analyse_maternal.sh
./run_analyse_paternal.sh
```

Each script calls `analyse_block.maternal/paternal.sh` for every chromosome, passing the chromosome name as an argument. You can also run a single chromosome directly:

```bash
bash analyse_block.maternal.sh chr1
```

Each variant passes through three layers of filtering:
- **PAN027 problematic regions** — regions where the mother's own assembly is unreliable
- **Grandparent problematic regions** — regions where the grandparent alt haplotype coordinate is unreliable (PAN010 for maternal, PAN011 for paternal)
- **Switch error-prone regions** — regions prone to haplotype switch errors in the pangenome

**Inheritance model:**
- **Maternal:** variant present in PAN010 (grandmother), absent in PAN028 (granddaughter) → candidate de novo in PAN027
- **Paternal:** variant present in PAN011 (grandfather), absent in PAN028 (granddaughter) → candidate de novo in PAN027

**Pipeline steps inside the template script:**

| Step | Description |
|------|-------------|
| 1 | Extract variants in de novo candidate blocks where grandparent carries the allele and granddaughter is REF (`bcftools` + AWK genotype check) |
| 2 | Build a BED of REF positions and run `vg find` on the GBZ graph (single call for all positions) |
| 3 | Map each REF position to the grandparent alt haplotype coordinate (`vg_pos_lookup.py`) |
| 4 | Insert `POS_PAN010/POS_PAN011` and `CONTIG_PAN010/CONTIG_PAN011` tags into the VCF INFO field |
| 5 | Trim unused ALT alleles (`bcftools view -a`) and remove exact duplicates (`bcftools norm -d exact`). Save unfiltered output. |
| 6 | Exclude variants in problematic/SW-prone regions (see below) |

**Filtering in STEP 6 — three exclusion criteria (OR logic):**

| Variable | What it excludes | BED file |
|----------|-----------------|----------|
| `EXCLUDE_A` | REF position overlaps PAN027 problematic regions | `problematic.PAN027.bed` |
| `EXCLUDE_B` | Grandparent alt position overlaps PAN010/PAN011 problematic regions | `problematic.PAN010.bed` / `problematic.PAN011.bed` |
| `EXCLUDE_C` | REF position overlaps switch error-prone regions | `sw_prone_regions.bed` |

Outputs per chromosome (in `results_final_filter_mother_gp/chr{N}_{mat|pat}/`):
- `*.wave.filtered.vcf` — after STEP 6 (final)
- `*.wave.unfiltered.vcf` — after STEP 5, before STEP 6

---

### Step 3 – Easy/difficult split (`filter_vcf_into_easy_and_difficult.sh`)

Splits each filtered VCF from Step 2 into two sets based on genomic region annotations:

- **easy** — variants intersecting `PAN027.v1.1.easy.bed`, with overlap removed (see below)
- **difficult** — variants intersecting `PAN027.v1.1.difficult.bed`

Variants that fall in both easy and difficult regions (overlap) are assigned to **difficult only** — they are subtracted from the easy set using `bcftools isec -C`. This ensures the easy and difficult sets are mutually exclusive.

Per-chromosome output (in `results_final_filter_mother_gp_easy_difficult/chr{N}_{mat|pat}/`):

| File | Description |
|------|-------------|
| `*.easy.raw.vcf.gz` | Easy variants before overlap removal |
| `*.easy.vcf.gz` | Easy variants after overlap removal (final) |
| `*.difficult.vcf.gz` | Difficult variants (includes overlap) |

The script skips chromosomes where all three output files already exist.

---

### Step 4 – Statistics and visualization (`run_make_stats_all.sh`)

Runs `make_stats.per_chr.py` and `create_chrom_map.py` for every chromosome and haplotype.

**`make_stats.per_chr.py`** — per-chromosome statistics:
- Classifies variants by `TYPE` and `LEN` from VCF INFO (with sequence-length fallback)
- MNP variants are decomposed: each base counted as one SNV
- Computes mutation counts and rates per bp (total block length and effective length)
- Effective block length = total block length minus problematic and SW-prone regions
- Outputs: `mutation_statistics.tsv`, `indel_lengths.tsv`, `heatmap_positions.bed`

**`create_chrom_map.py`** — per-chromosome visualization:
- Plots variant positions as vertical lines color-coded by type
- Analyzed blocks shown in grey
- Two styles: alpha overlay (all types on one track) and stacked tracks

---

## Variant Classification

| Category | Criterion |
|----------|-----------|
| SNV | `TYPE=snp`; or `TYPE=mnp` — each base counted as one SNV |
| Indel 1 bp | `TYPE=ins` or `del`, `LEN=1` |
| Indel 2–49 bp | `TYPE=ins` or `del`, `2 ≤ LEN ≤ 49` |
| SV | `TYPE=ins` or `del`, `LEN ≥ 50` |
| Tandem repeat | Overlap with TRF annotations (`PAN027_HiFi_element_polished.v1.0.unit1_6.trf.bed`) — additional label, may overlap any category |

When `TYPE=.` or `LEN=0`, values are derived from REF/ALT sequence lengths.

---

## Output Folder Structure

```
results_final_filter_mother_gp/
├── chr1_mat/
│   ├── PAN027.chr1.project_PAN027_mat.wave.filtered.vcf
│   ├── PAN027.chr1.project_PAN027_mat.wave.unfiltered.vcf
│   ├── mutation_statistics.tsv
│   ├── indel_lengths.tsv
│   ├── heatmap_positions.bed
│   ├── chromosome_map_alpha.png
│   └── chromosome_map_stacked.png
├── chr1_pat/
│   └── ...
└── ...

results_final_filter_mother_gp_easy_difficult/
├── chr1_mat/
│   ├── PAN027.chr1.project_PAN027_mat.wave.filtered.easy.raw.vcf.gz
│   ├── PAN027.chr1.project_PAN027_mat.wave.filtered.easy.vcf.gz
│   └── PAN027.chr1.project_PAN027_mat.wave.filtered.difficult.vcf.gz
└── ...
```

---

## Dependencies

- Cactus-pangenome v3.1.4 (via Podman)
- `vg` (for `vg find`, `vg paths`, `vg view`)
- `bcftools`
- `bedtools`
- `bgzip`, `tabix`
- Python 3 with: `pandas`, `matplotlib`
