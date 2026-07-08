#!/usr/bin/env python3
# Copyright (c) 2026, Dušan Slúka <xsluka@fi.muni.cz>
# SPDX-License-Identifier: BSD-3-Clause

"""
refine_pangenome_calls.py

Post-process the raw pangenome breakpoint calls (produced by
``parse_untangle.py``) with a bisection refinement that narrows each
call from ~100 kb–1 Mb down to ~1 kb.

Why this step is needed
-----------------------
``odgi untangle`` reports haplotype switches at the nearest structural
variation (SV) that lets the graph resolve which parental path the
child's haplotype follows. That anchor is often hundreds of kb away
from the true SNV-level identity flip. The refinement here re-anchors
each call to the SNV-level flip using minimap2 probe alignments against
a combined parental reference.

Inputs
------
- A pangenome switches.tsv (produced by ``parse_untangle.py``, e.g.
  ``PAN027_hap1.switches.tsv``).
- The child haplotype FASTA (matching the ``child-id`` in the switches
  TSV).
- A combined parent reference FASTA with contigs prefixed ``hap1|`` /
  ``hap2|`` (build once with ``samtools faidx`` + ``awk`` — see the
  ``config/COMMANDS.md`` document for the exact recipe).

Per chromosome, the script:
  1. Extracts the chr-only slice of the parent reference (samtools faidx)
  2. Extracts the chr-only slice of the child haplotype
  3. Synthesises candidate_calls.tsv rows (one per bp), expanded by ±N bp
  4. Invokes ``switch_coordinates_aproximation.py`` on the chr-restricted refs
  5. Concatenates refined coordinates into one output TSV
"""
import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path


HAP_TO_INTERNAL = {"hap1": "MATERNAL", "hap2": "PATERNAL"}


def slice_chr_fasta(src_fa: Path, contig_pattern: str, dst_fa: Path):
    """Use samtools faidx to write a single-contig FASTA matching pattern."""
    # Find matching contig name from .fai
    fai = src_fa.with_suffix(src_fa.suffix + ".fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(src_fa)], check=True)
    contigs = []
    with fai.open() as f:
        for line in f:
            name = line.split("\t", 1)[0]
            if re.search(contig_pattern, name):
                contigs.append(name)
    if not contigs:
        return []
    cmd = ["samtools", "faidx", str(src_fa)] + contigs
    with dst_fa.open("w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    subprocess.run(["samtools", "faidx", str(dst_fa)], check=True)
    return contigs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--switches", required=True, help="pangenome switches.tsv input")
    ap.add_argument("--child-fasta", required=True, help="full child haplotype FASTA")
    ap.add_argument("--mother-ref-fa", required=True,
                    help="combined parent FASTA (hap1|/hap2| prefixed)")
    ap.add_argument("--out", required=True, help="output refined TSV")
    ap.add_argument("--tmp-dir", default=None,
                    help="working directory (default: alongside output)")
    ap.add_argument("--expand", type=int, default=1_000_000,
                    help="expand each bp by ±N bp before bisection (default 1 Mb)")
    ap.add_argument("--iters", type=int, default=12)
    ap.add_argument("--probe-start", type=int, default=300_000)
    ap.add_argument("--probe-min", type=int, default=2_000)
    ap.add_argument("--resolution", type=int, default=1_000)
    ap.add_argument("--min-conf", type=float, default=0.2)
    ap.add_argument("--preset", default="asm5")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--strategy", choices=["bisection", "trisection"],
                    default="bisection")
    ap.add_argument("--bracket-edges", action="store_true")
    ap.add_argument("--bracket-max-shifts", type=int, default=5)
    ap.add_argument("--bracket-probe", type=int, default=50000)
    args = ap.parse_args()

    repo = Path(__file__).resolve().parent.parent
    bisect_script = repo / "scripts" / "switch_coordinates_aproximation.py"
    if not bisect_script.exists():
        sys.exit(f"[ERROR] missing {bisect_script}")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(args.tmp_dir) if args.tmp_dir else out_path.parent / "refine_tmp"
    tmp.mkdir(parents=True, exist_ok=True)

    # Group switches by chromosome
    by_chr = {}
    tag = None
    with open(args.switches) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            tag = r["tag"]
            by_chr.setdefault(r["chr"], []).append(r)
    if tag is None:
        sys.exit(f"[ERROR] no rows in {args.switches}")

    print(f"[refine] {args.switches}: {sum(len(v) for v in by_chr.values())} bps "
          f"on {len(by_chr)} chromosomes")

    # Output header
    header_written = False
    with out_path.open("w") as out:
        for chrom in sorted(by_chr.keys(), key=lambda c: (len(c), c)):
            rows = by_chr[chrom]
            print(f"\n[refine] === {chrom}: {len(rows)} bps ===")

            # Slice chr-restricted FASTAs
            chr_ref = tmp / f"{chrom}_parent_ref.fa"
            chr_child = tmp / f"{chrom}_child.fa"
            ref_contigs = slice_chr_fasta(Path(args.mother_ref_fa),
                                          rf"\.{chrom}\.", chr_ref)
            child_contigs = slice_chr_fasta(Path(args.child_fasta),
                                            rf"\.{chrom}\.", chr_child)
            if not ref_contigs or not child_contigs:
                print(f"[refine] skip {chrom}: missing contigs")
                continue

            # Synthesise candidate_calls.tsv
            cc = tmp / f"{chrom}_candidate_calls.tsv"
            with cc.open("w") as cf:
                cf.write("tag\tchr\tcandidate_bp\tleft_probe\tright_probe\t"
                         "left_winner\tright_winner\tcall\n")
                for r in rows:
                    bp = int(r["bp_est_0based"])
                    lo = max(0, bp - args.expand)
                    hi = bp + args.expand
                    lw = HAP_TO_INTERNAL.get(r["left_hap"])
                    rw = HAP_TO_INTERNAL.get(r["right_hap"])
                    if lw is None or rw is None:
                        print(f"[refine] skip bp {bp}: unknown hap "
                              f"{r['left_hap']}/{r['right_hap']}")
                        continue
                    cf.write(f"{r['tag']}\t{chrom}\t{lo}-{hi}\t"
                             f"dummy_left\tdummy_right\t{lw}\t{rw}\tSWITCH\n")

            # Run bisection
            chr_tmp = tmp / f"{chrom}_bisect"
            chr_tmp.mkdir(exist_ok=True)
            chr_out = tmp / f"{chrom}_refined.tsv"

            cmd = [
                sys.executable, str(bisect_script),
                "--candidate-calls", str(cc),
                "--child-fasta", str(chr_child),
                "--mother-ref-fa", str(chr_ref),
                "--hap1-pattern", "hap1|",
                "--hap2-pattern", "hap2|",
                "--iters", str(args.iters),
                "--probe-start", str(args.probe_start),
                "--probe-min", str(args.probe_min),
                "--resolution", str(args.resolution),
                "--min-conf", str(args.min_conf),
                "--preset", args.preset,
                "--threads", str(args.threads),
                "--secondary", "no",
                "--cigar",
                "--strategy", args.strategy,
                "--tmp-dir", str(chr_tmp),
                "--out", str(chr_out),
            ]
            if args.bracket_edges:
                cmd += ["--bracket-edges",
                        "--bracket-max-shifts", str(args.bracket_max_shifts),
                        "--bracket-probe", str(args.bracket_probe)]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                print(f"[refine] bisection failed on {chrom} (exit "
                      f"{proc.returncode})\nstderr:\n{proc.stderr[-1000:]}")
                continue

            # Read refined output and append to global out
            with chr_out.open() as rf:
                hdr = rf.readline().rstrip("\n")
                if not header_written:
                    out.write(hdr + "\n")
                    header_written = True
                for line in rf:
                    out.write(line)

    print(f"\n[OK] {out_path}")


if __name__ == "__main__":
    main()
