# Pangenome-based meiotic recombination breakpoint pipeline

A pipeline for detecting meiotic crossover breakpoints in haplotype-resolved
trio assemblies using a pangenome graph representation. Developed as part of
the PAN027 recombination analysis for the washu-pedigree dataset.

- **Input:** two phased haplotype FASTAs each for child, mother, and father
  (six diploid haplotypes total, PanSN-formatted headers).
- **Output:** breakpoints where the child's haplotype path switches between
  the two parental haplotypes, with bp-level coordinates in TSV, BED, and PAF.

## Method

The pipeline has three stages:

1. **Pangenome graph construction** — `pggb` builds a single variation graph
   per chromosome over the six trio haplotypes (`wfmash` → `seqwish` →
   `smoothxg`).
2. **Untangling** — `odgi untangle` walks each child path through the graph
   and reports, for every short segment, which parental path shares the most
   graph nodes. Contiguous segments assigned to the same parental haplotype
   are merged; each haplotype flip is a candidate crossover.
3. **Local refinement** (this work, `refine_pangenome_calls.py`) — narrows
   each candidate to ~1 kb precision with minimap2-based probe alignment. A
   pre-phase (*edge bracketing*) confirms that the switch is actually
   contained in the search window and drops candidates that turn out to be
   false positives from repetitive-region graph ambiguity.

## Repository layout

```
pangenomes/
├── README.md               # this file
├── breakpoints/            # final PAN027 breakpoint calls
│   ├── PAN027_maternal.breakpoints.tsv
│   └── PAN027_paternal.breakpoints.tsv
├── config/
│   ├── PAN027_final_params.json    # exact pipeline config used
│   ├── PARAMS_SUMMARY.md           # concise one-file parameter summary
│   └── COMMANDS.md                 # step-by-step CLI to reproduce
├── scripts/
│   ├── parse_untangle.py                  # odgi BED → switches.tsv
│   ├── refine_pangenome_calls.py          # bracket + bisection refinement
│   ├── switch_coordinates_aproximation.py # bracket + bisection core
│   └── plot_pangenome_karyogram.py        # karyogram of a run
└── figures/                # example outputs


```

## Install

```bash
conda install -c bioconda -c conda-forge pggb=0.7.4 minimap2 samtools matplotlib
# odgi: build from source (>= 0.9.4), bioconda binaries crash inside `untangle`
git clone --recursive https://github.com/pangenome/odgi.git
cd odgi && cmake -H. -Bbuild && cmake --build build -- -j 16
# add <odgi>/bin to $PATH before running the untangle step
```

Full parameter dump used for the PAN027 v1.0 run is in
[`config/PAN027_final_params.json`](config/PAN027_final_params.json).

## Final parameters (PAN027 v1.0 result)

### pggb (per-chromosome)

| Parameter | Value |
| --- | ---: |
| `-s` (wfmash segment length) | 50 000 bp |
| `-p` (wfmash percent identity) | 98 |
| `-k` (seqwish min match) | 23 |
| `-n` (haplotype count) | 6 |
| `-t` (threads) | 16 |

### odgi untangle

| Parameter | Value |
| --- | ---: |
| `--merge-dist` | 50 000 bp |
| `-n` (n-best) | 1 |
| Jaccard identity threshold (parser) | 0.90 |

Built from `odgi` source at commit tag ≥ 0.9.4; the bioconda binaries at
0.9.2 / 0.9.4 crash during `untangle` for large graphs.

### Refinement (bracket + bisection, minimap2 probes)

| Parameter | Value |
| --- | ---: |
| Window around each pangenome bp | ± 1 000 000 bp |
| Bisection iters | 12 |
| Probe start length | 300 000 bp |
| Probe min length | 2 000 bp |
| Convergence resolution | 1 000 bp |
| Min probe-alignment confidence | 0.20 |
| Bracket max window shifts | 5 |
| Bracket edge-probe length | 50 000 bp |
| Min clamped probe length | 500 bp |
| minimap2 preset | `asm5`, `--secondary=no`, `-c` |

The refinement wall-clock is ≈ 45 min on a shared 4-CPU node with the
mother and father branches running in parallel.

## Reproducing the final run

Full step-by-step commands are in [`config/COMMANDS.md`](config/COMMANDS.md).
Short version, assuming the six PAN027-trio haplotype FASTAs and a
combined parental reference `mother_ref.fa` (with `hap1|` / `hap2|`
prefixed contigs):

```bash
# 1. Build the pangenome graph (per-chromosome)
pggb -i chr1.fa -o graph/chr1 \
     -s 50000 -p 98 -k 23 -n 6 -t 16

# 2. Untangle: PAN027 hap1 (child mother-derived) vs PAN010 (mother)
odgi untangle -i graph/chr1/chr1.final.og \
     -R target_paths.txt -Q query_paths.txt \
     -t 4 --merge-dist 50000 > untangle/PAN027_hap1_vs_PAN010.bed

# 3. Parse untangle ut into switches
python3 scripts/parse_untangle.py \
     --bed  untangle/PAN027_hap1_vs_PAN010.bed \
     --child-id PAN027 --parent-id PAN010 --hap hap1 \
     --min-identity 0.90 \
     --out  untangle/PAN027_hap1.switches.tsv

# 4. Refine coordinates with bracket + bisection
python3 scripts/refine_pangenome_calls.py \
     --switches       untangle/PAN027_hap1.switches.tsv \
     --child-fasta    PAN027_hap1.fasta \
     --mother-ref-fa  mother_ref.fa \
     --iters 12 --probe-min 2000 --resolution 1000 \
     --bracket-edges --bracket-max-shifts 5 --bracket-probe 50000 \
     --threads 4 \
     --out            refined/PAN027_hap1.refined.tsv

# 5. Plot the karyogram of the resulting BED
python3 scripts/plot_pangenome_karyogram.py \
     --run-dir refined/ --fai PAN027_hap1.fasta.fai \
     --out figures/pangenome_final.png
```

Repeat steps 2–4 for the paternal branch (`hap2` vs the father).

## Output format

`breakpoints/PAN027_{maternal,paternal}.breakpoints.tsv` columns:

| Column | Meaning |
| --- | --- |
| `chr` | `chr1`, `chr2`, … |
| `bp_est_1based` | best point estimate of the crossover (1-based) |
| `bp_lo`, `bp_hi` | uncertainty interval (bisection convergence bounds) |
| `conf` | probe-agreement confidence at convergence |
| `left_hap`, `right_hap` | flanking parental haplotypes (`hap1` / `hap2`) |
| `status` | see below |


Rows with `status ∈ {LIKELY_FP, EMPTY_PROBE}` are dropped for downstream
analysis.

### `status` values

| Status | Meaning | Use |
| --- | --- | --- |
| `OK` | Bisection converged to a ≤ 1 kb interval — a well-localized breakpoint. | **valid** |
| `ACTIVE` | A real breakpoint whose interval did not fully close to 1 kb within the iteration budget (e.g. sparse discriminating variants), but the switch is confirmed and localized to a few kb. | **valid** |
| `LIKELY_FP` | Edge bracketing never found a HAP1↔HAP2 transition anywhere in the search span — no crossover exists in the window, so the raw call was a false positive (graph ambiguity in a repeat). | dropped |
| `EMPTY_PROBE` | The candidate sits so close to a chromosome end that the clamped probe collapses below the minimum length and cannot be aligned — the breakpoint cannot be validated (see Limitations). | dropped |

Valid breakpoints are those with `status ∈ {OK, ACTIVE}`.


## Limitations

- **Breakpoints very close to a telomere cannot be validated.** The refinement
  places a probe on each side of the candidate to check which parental
  haplotype it matches. Within roughly a probe-minimum length (~500 bp to a
  few kb) of a chromosome end there is not enough sequence on the outer side to
  extract a probe, so such candidates are flagged `EMPTY_PROBE` and dropped.
  Because paternal crossovers are biologically enriched near telomeres, this can
  under-represent the most telomere-proximal paternal breakpoints.
- **Chromosome-scoped graphs.** Graphs are built per chromosome, so chrX↔chrY
  pseudoautosomal recombination is not detectable in the default run; a combined
  chrX+chrY graph is needed for the PAR.
- **Doublet artifacts.** The graph occasionally emits two calls a few kb apart
  (a switch and an immediate switch-back); consecutive calls < ~200 kb apart on
  the same haplotype should be treated as a single crossover.

