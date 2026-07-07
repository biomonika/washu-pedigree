# Final parameters summary — PAN027 v1.0

Concise one-file summary of parameters used for the final published run.
For a full JSON dump see `PAN027_final_params.json`, and for step-by-step
CLI commands see `COMMANDS.md`.

## Pangenome graph (pggb)

- `--segment-length 50000`
- `--percent-identity 98`
- `--kmer 23`
- `--n-haplotypes 6`
- Per-chromosome graphs (not whole-genome).

## Untangling (odgi untangle)

- `--merge-dist 50000`
- Best hit only (`--n-best 1`, default).
- Post-parse Jaccard filter: >= 0.90.

## Refinement (bracket + bisection)

- Initial search window: ±1 Mb around each raw pangenome bp
- Bisection: 12 iters, probe shrinks 300 kb → 2 kb, stop at 1 kb window
- Edge bracketing: 50 kb edge probes, up to 5 window shifts
- Window clamped to chromosome bounds (probes may be asymmetric near
  telomeres, minimum probe length 500 bp).
- minimap2 preset `asm5`, `--secondary=no`, `-c`.
