# Commands to reproduce the final PAN027 pangenome run

The pangenome graph is built per-chromosome, then untangled twice per
chromosome (child hap1 vs mother, child hap2 vs father), then refined.

## 0. Prerequisites

- Six haplotype FASTAs in a common directory:

```
data/
  PAN010_hap1.fasta   # grandmother
  PAN010_hap2.fasta
  PAN011_hap1.fasta   # grandfather
  PAN011_hap2.fasta
  PAN027_hap1.fasta   # mother — child of the (PAN010, PAN011) meiosis
  PAN027_hap2.fasta
```

- `pggb ≥ 0.7.4`, `odgi` **built from source at ≥ 0.9.4** (bioconda binaries
  crash during `untangle` on large graphs), `minimap2`, `samtools`,
  Python 3.10 with `matplotlib` (install commands are in the README).

- A combined parental reference for the refinement, with `hap1|` and
  `hap2|` prefixed contigs. For the mother side, build it once:

```bash
awk '/^>/{print ">hap1|"substr($0,2); next} {print}' PAN010_hap1.fasta \
    > mother_ref.fa
awk '/^>/{print ">hap2|"substr($0,2); next} {print}' PAN010_hap2.fasta \
    >> mother_ref.fa
samtools faidx mother_ref.fa
```

Repeat with PAN011 for `father_ref.fa`.

## 1. Combined PanSN input (per chromosome)

Split each haplotype FASTA by chromosome and rename headers to PanSN
(`SAMPLE#HAP#CHR`), then concatenate all six haplotypes' chromosome-N
sequences into `chrN.fa`. Example for chr1:

```bash
for s in PAN010 PAN011 PAN027 ; do
  for h in 1 2 ; do
    awk -v s=$s -v h=$h '/^>/{print ">"s"#"h"#chr1"; next} {print}' \
        data/${s}_hap${h}.chr1.fasta
  done
done > per_chr/chr1.fa
samtools faidx per_chr/chr1.fa
```

Do the same for chr2 … chr22 (and chrX if desired).

## 2. Build the pangenome graph

Per chromosome:

```bash
pggb \
  -i per_chr/chr1.fa \
  -o graph/chr1 \
  -s 50000 \
  -p 98 \
  -k 23 \
  -n 6 \
  -t 16
```

The output graph is `graph/chr1/chr1.fa.<hashes>.smooth.final.og`.

## 3. Untangle each parent branch

Mother side (child hap1 vs both mother haplotypes):

```bash
echo "PAN010#1#chr1" >  target_mother.txt
echo "PAN010#2#chr1" >> target_mother.txt
echo "PAN027#1#chr1" >  query_hap1.txt

odgi untangle \
  -i graph/chr1/chr1.*.smooth.final.og \
  -R target_mother.txt \
  -Q query_hap1.txt \
  -t 4 \
  --merge-dist 50000 \
  > untangle/PAN027_hap1_vs_PAN010.chr1.bed
```

Do the same for the father side (`hap2` vs both PAN011 haplotypes) and
concatenate the per-chromosome BEDs into one file each.

## 4. Parse into switches

```bash
python3 scripts/parse_untangle.py \
  --bed          untangle/PAN027_hap1_vs_PAN010.bed \
  --child-id     PAN027 \
  --parent-id    PAN010 \
  --hap          hap1 \
  --min-identity 0.90 \
  --out          untangle/PAN027_hap1.switches.tsv
```

Repeat with `--hap hap2 --parent-id PAN011` for the paternal side.

## 5. Refine coordinates

```bash
python3 scripts/refine_pangenome_calls.py \
  --switches       untangle/PAN027_hap1.switches.tsv \
  --child-fasta    data/PAN027_hap1.fasta \
  --mother-ref-fa  mother_ref.fa \
  --out            refined/PAN027_hap1.refined.tsv \
  --tmp-dir        refined/tmp_mother \
  --expand         1000000 \
  --iters          12 \
  --probe-start    300000 \
  --probe-min      2000 \
  --resolution     1000 \
  --min-conf       0.20 \
  --bracket-edges \
  --bracket-max-shifts 5 \
  --bracket-probe  50000 \
  --preset         asm5 \
  --threads        4
```

Repeat with `--mother-ref-fa father_ref.fa`,
`--child-fasta data/PAN027_hap2.fasta`,
`--out refined/PAN027_hap2.refined.tsv` for the paternal side.

Refined coordinates keep rows tagged `OK` and `ACTIVE`; rows tagged
`LIKELY_FP` (no HAP1↔HAP2 transition observable in the search span) or
`EMPTY_PROBE` (probe collapsed below 500 bp near a chromosome end) should
be dropped for downstream analysis.

## 6. Plot

```bash
# Optional: combine both branches into a single BED for plotting
{ tail -n +2 refined/PAN027_hap1.refined.tsv ; \
  tail -n +2 refined/PAN027_hap2.refined.tsv ; } | \
awk -F'\t' 'BEGIN{OFS="\t"}
  $NF!="LIKELY_FP" && $NF!="EMPTY_PROBE" {
    print $2, $5, $6, $1"_SWITCH", int($8*1000), "+"
  }' > refined/crossovers.bed

python3 scripts/plot_pangenome_karyogram.py \
  --run-dir refined/ \
  --fai     data/PAN027_hap1.fasta.fai \
  --out     figures/pangenome_final.png
```

The plotter auto-detects `crossovers.bed` under the run directory.
