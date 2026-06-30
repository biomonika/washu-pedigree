# Scripts for using centrodip to call CDRs

Tool versions: `minimap2 2.30-r1287`, `modkit v0.5.0`, `centrodip v0.1.0-pre1` (bundled as `centrodip-0.1.0-pre1.tar.gz`). Throughout, reads are filtered to primary alignments with MAPQ ≥ 10, and the active-alpha array is the `active_hor` cenSat annotation.

## `centrodip_workflow.sh`
Calls CDRs on each sample independently (no cross-sample realignment). Holds in-script paths to the cenSat BEDs, ONT methylation BAMs, and polished reference FASTAs.

**For each sample's set of BAM, cenSat BED, and reference FASTA:**
1. Realign reads to that sample's own matched polished reference with `minimap2 -ax map-ont -I 16G --eqx --cs -Y -L -y -p 0.5 -k 17` (keep MAPQ ≥ 10).

2. Aggregate mCpG with `modkit pileup --cpg --combine-strands --mod-thresholds m:0.8 --filter-threshold C:0.5`, then keep only positions with coverage ≥ 5.

3. Subset the cenSat to active-alpha (`active_hor`) and run `centrodip` on the filtered modkit pileup to call CDRs.

---

# Realignment analysis and visualization

## `pedigree_centromere_realignment.map_ont.sh` & `pedigree_centromere_realignment.lr_hq.sh`
Holds in-script paths to cenSat annotations, BAMs, and polished reference FASTAs. The two scripts are identical except for the `minimap2` preset and the minimum read length:
- `map_ont` — `-ax map-ont`, minimum read length 10 kb
- `lr_hq` — `-ax lr:hq`, minimum read length 100 kb

**For each sample's set of BAM, cenSat BED, and reference FASTA:**

1. Align reads to the matched assembly.

2. For each chromosome:
   - Loop through haplotypes:
     - Extract the active-alpha reads for that chromosome/haplotype combination
     - Realign those reads to PAN027 (diploid) with `minimap2` (same flags as above)

3. For each PAN027 haplotype of that chromosome:
   - Evaluate mapping efficiency (fraction of MAPQ ≥ 10 reads retained after realignment to the PAN027 haplotype)
   - Aggregate methylation with `modkit pileup --cpg --combine-strands --mod-thresholds m:0.8 --filter-threshold C:0.5`
   - Call CDRs with `centrodip`

## `visualization.13101294.ipynb`
Notebook for visualizing inheritance patterns in centromeric methylation, local identity, HORhaps (k=9), and de novo candidate variants. (13101294 is the `SLURM_JOB_ID` of the `map_ont` realignment run.)
