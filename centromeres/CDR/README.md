# Scripts for using centrodip to call CDRs

## `centrodip_workflow.sh`
Workflow which contains in paths of cenSat, BAMs, and polished reference FASTAs.  

**For each samples set of BAMs, cenSat BED, and reference FASTA:**
1. Realigns reads to matched polished reference genome

2. call `modkit` to aggregate mCpG
    * subset to just active-alpha

3. Run `centrodip` on modkit output and active-alpha cenSat

---

# Realignment analysis and visualization

## `pedigree_centromere_realignment.map_ont.sh` & `pedigree_centromere_realignment.lr_hq.sh`
Contains paths to cenSat annotations, BAMs, and polished reference FASTAs

**For each samples set of BAMs, cenSat BED, and reference FASTA:**

1. Perform alignment to the matched assembly.

2. For each chromosome:
   - Loop through haplotypes:
     - Extract reads for the chromosome/haplotype combination
     - Realign those reads to PAN027 (diploid)

3. For each PAN027 haplotype of that chromosome:
   - Evaluate mapping efficiency
   - Call aggregated methylation
   - Call CDRs

## `visualization.13101294.ipynb`
Notebook for visualizing inheritance patterns in centromeric methylation, local identity, HORhaps, and mutations. (13101294 is the SLURM_JOB_ID from realignment)