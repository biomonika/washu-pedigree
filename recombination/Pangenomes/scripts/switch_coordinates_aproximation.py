#!/usr/bin/env python3
# Copyright (c) 2026, Dušan Slúka <xsluka@fi.muni.cz>
# SPDX-License-Identifier: BSD-3-Clause

"""
switch_coordinates_aproximation.py

Narrows down the position of each candidate crossover to ~1 kb resolution
via minimap2-based probe alignment.

Called by ``refine_pangenome_calls.py`` — a raw pangenome switch is
expanded to a ±1 Mb search window, this script confirms the switch is
inside that window (edge bracketing) and then bisects it. Rows that
cannot be bracketed are flagged ``LIKELY_FP``; rows whose probes fall
past a chromosome end and cannot be clamped to the minimum size are
flagged ``EMPTY_PROBE``. Both are dropped by downstream filtering.

Inputs
------
- ``--candidate-calls``   TSV of pangenome bp candidates (produced by
                          ``refine_pangenome_calls.py``). Columns
                          ``tag chr candidate_bp left_winner right_winner
                          call``; each ``candidate_bp`` is a
                          ``LO-HI`` window.
- ``--child-fasta``       Full child haplotype FASTA (must include the
                          chromosomes referenced in ``candidate-calls``).
- ``--mother-ref-fa``     Combined parent reference (``hap1|`` /
                          ``hap2|`` prefixed contigs).

Edge bracketing (pre-phase)
---------------------------
Before bisecting, two ``--bracket-probe`` bp probes at the interval
edges verify that the true switch is bracketed. If both edges align to
the same haplotype, the window is walked in that direction; if walking
would exit the chromosome, the window is clamped and either the walk
continues on asymmetric probes or the candidate is flagged
``LIKELY_FP`` (no HAP1↔HAP2 transition anywhere in the reachable span).

Bisection (main loop)
---------------------
For each candidate, place two probes of length *P* on either side of
the midpoint ``m = (lo + hi) / 2``:

  1. Both align to ``left_hap``  → ``lo = m`` (switch is in the right half)
  2. Both align to ``right_hap`` → ``hi = m`` (switch is in the left half)
  3. Different flanks            → ``m`` is the boundary; tighten to
                                    ``[m − P/2, m + P/2]``

The probe length is halved with each iteration
(``P = max(probe_min, probe_start / 2^(i-1))``) so probes never cross
the boundary they are trying to localise. Iteration stops when the
interval is ≤ ``--resolution`` or after ``--iters`` rounds. Optional
``--strategy trisection`` places two test points at 1/3 and 2/3 of the
interval for potentially faster convergence per iteration.

The classical bisection method (Burden, Faires & Burden,
*Numerical Analysis*, 10th ed., 2016, ch. 2.1) is adapted here to
discrete probe evidence rather than a continuous function root.
"""
import argparse
import csv
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, List, Optional


def which(exe: str) -> Optional[str]:
    return shutil.which(exe)


def find_contig_from_fai(fai_path: str, chr_name: str) -> str:
    """Resolve a plain chromosome name (``chr1``, ``chr22`` …) to the
    full FASTA contig name in ``fai_path`` (e.g. ``PAN027.chr1.maternal``)."""
    c = chr_name.replace("chr", "")
    pat = re.compile(rf"\.chr{re.escape(c)}\.", re.IGNORECASE)
    with open(fai_path, "r") as f:
        for line in f:
            contig = line.split("\t", 1)[0]
            if pat.search(contig):
                return contig
    return ""


def contig_length(fai_path: str, contig: str) -> int:
    """Return the length of a contig from the FAI (second column). 0 if not found."""
    with open(fai_path, "r") as f:
        for line in f:
            parts = line.split("\t")
            if parts and parts[0] == contig:
                try:
                    return int(parts[1])
                except (ValueError, IndexError):
                    return 0
    return 0


def faidx_region(fasta: str, region: str) -> str:
    # samtools faidx uses 1-based inclusive coordinates; a start of 0 is
    # rejected as "Zero length sequence" and returns empty. Rewrite any
    # `contig:0-END` to `contig:1-END` so telomere-adjacent probes work.
    if ":0-" in region:
        head, tail = region.split(":0-", 1)
        region = f"{head}:1-{tail}"
    cmd = ["samtools", "faidx", fasta, region]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"samtools faidx failed: {' '.join(cmd)}\n{p.stderr}")
    return p.stdout


def fasta_to_seq(fasta_text: str) -> str:
    out = []
    for ln in fasta_text.splitlines():
        if ln.startswith(">"):
            continue
        ln = ln.strip()
        if ln:
            out.append(ln)
    return "".join(out)


def write_fasta_record(fh, name: str, seq: str, width: int = 60):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i+width] + "\n")


def score_paf_by_patterns(paf_path: Path, hap1_pat: str, hap2_pat: str) -> Dict[str, Tuple[int, int, str, float]]:
    """Scores a minimap2 PAF the same way ``score_probes_paf.py`` does —
    sums alignment lengths against hap1 and hap2 of the parent and
    chooses the larger as the winner.

    The confidence is the normalised difference
    ``|hap1 − hap2| / (hap1 + hap2)``; values close to 0 mean the
    probe could not discriminate and bisection should stop.

    Returns a mapping ``probe_id -> (hap1_aln, hap2_aln, winner, conf)``
    where ``winner`` is one of ``MATERNAL``, ``PATERNAL``, ``TIE`` or
    ``MISSING``."""
    a: Dict[str, int] = {}
    b: Dict[str, int] = {}

    with paf_path.open("r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            q = parts[0]
            t = parts[5]
            aln = int(parts[10])
            if hap1_pat in t:
                a[q] = a.get(q, 0) + aln
            elif hap2_pat in t:
                b[q] = b.get(q, 0) + aln

    out: Dict[str, Tuple[int, int, str, float]] = {}
    keys = set(a) | set(b)
    for k in keys:
        x = a.get(k, 0)
        y = b.get(k, 0)
        tot = x + y
        if tot == 0:
            out[k] = (0, 0, "MISSING", 0.0)
            continue
        if x == y:
            out[k] = (x, y, "TIE", 0.0)
            continue
        winner = "MATERNAL" if x > y else "PATERNAL"
        conf = abs(x - y) / tot
        out[k] = (x, y, winner, conf)
    return out


@dataclass
class Cand:
    tag: str
    chr: str
    lo: int
    hi: int
    left_exp: str
    right_exp: str
    call: str

    # tracking
    iters_done: int = 0
    last_mid: int = -1
    last_test_points: list = None
    last_left: str = "NA"
    last_right: str = "NA"
    last_conf: float = 0.0
    status: str = "ACTIVE"


def parse_candidate_calls(path: str, include_non_switch: bool) -> List[Cand]:
    out: List[Cand] = []
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            tag = row["tag"]
            chr_ = row["chr"]
            call = row["call"]
            if (call != "SWITCH") and (not include_non_switch):
                continue

            s, e = row["candidate_bp"].split("-", 1)
            lo = int(float(s))
            hi = int(float(e))

            left_exp = row["left_winner"]
            right_exp = row["right_winner"]

            out.append(Cand(tag=tag, chr=chr_, lo=lo, hi=hi,
                            left_exp=left_exp, right_exp=right_exp, call=call))
    return out


def run_minimap(mother_ref: str, probes_fa: Path, out_paf: Path,
               preset: str, threads: int, secondary: str, cigar: bool):
    cmd = ["minimap2", "-x", preset, "-t", str(threads), f"--secondary={secondary}"]
    if cigar:
        cmd.append("-c")
    cmd += [mother_ref, str(probes_fa)]
    with out_paf.open("w") as fo:
        p = subprocess.run(cmd, stdout=fo)
    if p.returncode != 0:
        raise SystemExit(f"[ERROR] minimap2 failed (exit {p.returncode})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate-calls", required=True, help="*.candidate_calls.tsv")
    ap.add_argument("--child-fasta", required=True)
    ap.add_argument("--mother-ref-fa", required=True)
    ap.add_argument("--hap1-pattern", default="maternal")
    ap.add_argument("--hap2-pattern", default="paternal")

    ap.add_argument("--iters", type=int, default=8)
    ap.add_argument("--probe-start", type=int, default=300000)
    ap.add_argument("--probe-min", type=int, default=2000)
    ap.add_argument("--resolution", type=int, default=1000)
    ap.add_argument("--min-conf", type=float, default=0.2)
    ap.add_argument("--include-non-switch", action="store_true")

    ap.add_argument("--preset", default="asm20")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--secondary", default="no")
    ap.add_argument("--cigar", action="store_true")

    ap.add_argument("--tmp-dir", required=True)
    ap.add_argument("--out", required=True)

    ap.add_argument("--strategy", choices=["bisection", "trisection"],
                    default="bisection",
                    help="Bisection halves the interval per iter (1 test point, 2 probes). "
                         "Trisection uses 2 test points (4 probes) at 1/3 and 2/3 of the "
                         "interval and can shrink to 1/3 per iter in the best case.")

    ap.add_argument("--bracket-edges", action="store_true",
                    help="Before bisection: probe interval edges to confirm the switch "
                         "is bracketed by the [lo, hi] window. If both edges match the "
                         "same haplotype, shift the window in that direction up to "
                         "--bracket-max-shifts times. After that, mark as LIKELY_FP "
                         "(probable pangenome false positive — no actual switch nearby).")
    ap.add_argument("--bracket-max-shifts", type=int, default=5,
                    help="Max window shifts during edge bracketing (default 5).")
    ap.add_argument("--bracket-probe", type=int, default=50000,
                    help="Probe length for edge tests (default 50 kb).")

    args = ap.parse_args()

    if not which("samtools"):
        raise SystemExit("[ERROR] samtools not in PATH")
    if not which("minimap2"):
        raise SystemExit("[ERROR] minimap2 not in PATH")

    if not Path(args.mother_ref_fa).exists():
        raise SystemExit(f"[ERROR] mother-ref-fa missing: {args.mother_ref_fa}")

    # Ensure FASTA index exists
    fai = args.child_fasta + ".fai"
    if not os.path.exists(fai):
        p = subprocess.run(["samtools", "faidx", args.child_fasta])
        if p.returncode != 0:
            raise SystemExit("[ERROR] samtools faidx failed for child-fasta")

    cands = parse_candidate_calls(args.candidate_calls, args.include_non_switch)
    if not cands:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["tag","chr","bp_est_0based","bp_est_1based","bp_lo","bp_hi",
                        "iters","conf","left_exp","right_exp","last_left","last_right","status"])
        print("[OK] no candidates selected")
        return

    tmp = Path(args.tmp_dir)
    tmp.mkdir(parents=True, exist_ok=True)

    chr2contig: Dict[str, str] = {}
    chr2length: Dict[str, int] = {}
    for c in cands:
        if c.chr not in chr2contig:
            contig = find_contig_from_fai(fai, c.chr)
            chr2contig[c.chr] = contig
            chr2length[c.chr] = contig_length(fai, contig) if contig else 0

    # Clamp obviously bad intervals to chromosome bounds
    for c in cands:
        if c.lo < 0:
            c.lo = 0
        chr_len = chr2length.get(c.chr, 0)
        if chr_len > 0 and c.hi > chr_len:
            c.hi = chr_len
        if c.hi <= c.lo:
            c.status = "BAD_INTERVAL"

    # Minimum probe length below which an edge probe is considered too short
    # to be informative (used near chromosome ends). Below this the edge
    # cannot vote — but LIKELY_FP is still decided per --bracket-max-shifts.
    MIN_PROBE_LEN = 500

    active = [c for c in cands if c.status == "ACTIVE"]

    # ── Edge bracketing phase ─────────────────────────────────────────────
    # Probes at [lo, lo+P] and [hi-P, hi] confirm whether the switch lies
    # between them. If both edges show the same haplotype, shift the window
    # in that direction (window-walking). After --bracket-max-shifts failed
    # attempts the candidate is marked LIKELY_FP — a pangenome call with no
    # actual haplotype switch nearby.
    if args.bracket_edges and active:
        bp = args.bracket_probe
        for shift_it in range(args.bracket_max_shifts + 1):
            if not active:
                break
            probes_fa = tmp / f"bracket{shift_it:02d}.probes.fa"
            paf = tmp / f"bracket{shift_it:02d}.paf"

            def _edge_regions(c):
                """Return ((L1, L2), (R1, R2)) clamped to [0, chr_len].
                Each edge probe is up to `bp` bp long; near chromosome ends
                it is silently shortened so refinement can still proceed."""
                chr_len = chr2length.get(c.chr, 0) or c.hi
                L1 = max(0, c.lo)
                L2 = min(chr_len, L1 + bp)
                R2 = min(chr_len, c.hi)
                R1 = max(0, R2 - bp)
                return (L1, L2), (R1, R2)

            with probes_fa.open("w") as outfa:
                for c in active:
                    contig = chr2contig.get(c.chr, "")
                    if not contig:
                        c.status = "NO_CONTIG"
                        continue
                    (L1, L2), (R1, R2) = _edge_regions(c)
                    # Only the (rare) case where a probe collapses to < 500 bp
                    # is treated as empty; otherwise we allow asymmetric probes.
                    if (L2 - L1) < MIN_PROBE_LEN or (R2 - R1) < MIN_PROBE_LEN:
                        c.status = "EMPTY_PROBE"
                        continue
                    left_seq = fasta_to_seq(faidx_region(args.child_fasta, f"{contig}:{L1}-{L2}"))
                    right_seq = fasta_to_seq(faidx_region(args.child_fasta, f"{contig}:{R1}-{R2}"))
                    if not left_seq or not right_seq:
                        c.status = "EMPTY_PROBE"
                        continue
                    lid = f"{c.tag}|{c.chr}|BRACKET|sh{shift_it}|L|{L1}-{L2}"
                    rid = f"{c.tag}|{c.chr}|BRACKET|sh{shift_it}|R|{R1}-{R2}"
                    write_fasta_record(outfa, lid, left_seq)
                    write_fasta_record(outfa, rid, right_seq)

            run_minimap(args.mother_ref_fa, probes_fa, paf,
                       preset=args.preset, threads=args.threads,
                       secondary=args.secondary, cigar=args.cigar)
            scores = score_paf_by_patterns(paf, args.hap1_pattern, args.hap2_pattern)

            next_active = []
            for c in active:
                if c.status != "ACTIVE":
                    continue
                (L1, L2), (R1, R2) = _edge_regions(c)
                lid = f"{c.tag}|{c.chr}|BRACKET|sh{shift_it}|L|{L1}-{L2}"
                rid = f"{c.tag}|{c.chr}|BRACKET|sh{shift_it}|R|{R1}-{R2}"
                lw = scores.get(lid, (0,0,"MISSING",0.0))[2]
                rw = scores.get(rid, (0,0,"MISSING",0.0))[2]

                if lw == c.left_exp and rw == c.right_exp:
                    # Correctly bracketed → leave window as-is, proceed to bisection
                    continue
                if lw == c.right_exp and rw == c.left_exp:
                    # Inverted bracketing (only valid if user swapped expected sides)
                    continue
                # Both edges show the same haplotype → shift window in that direction
                W = c.hi - c.lo
                chr_len = chr2length.get(c.chr, 0)
                if lw == c.left_exp and rw == c.left_exp:
                    new_lo = c.hi
                    new_hi = c.hi + W
                    if chr_len > 0:
                        new_hi = min(new_hi, chr_len)
                        new_lo = min(new_lo, max(0, new_hi - 1))
                    if new_hi <= new_lo:
                        # Cannot shift further right without exiting chromosome
                        c.status = "LIKELY_FP"
                        continue
                    c.lo, c.hi = new_lo, new_hi
                    next_active.append(c)
                elif lw == c.right_exp and rw == c.right_exp:
                    new_hi = c.lo
                    new_lo = max(0, c.lo - W)
                    if new_hi <= new_lo:
                        c.status = "LIKELY_FP"
                        continue
                    c.lo, c.hi = new_lo, new_hi
                    next_active.append(c)
                else:
                    # MISSING / TIE / inconsistent — try once more, else give up
                    if shift_it < args.bracket_max_shifts:
                        next_active.append(c)
                    else:
                        c.status = "LIKELY_FP"
            active = next_active

        # Any candidate still in 'active' after max shifts failed to bracket
        for c in active:
            c.status = "LIKELY_FP"
        active = [c for c in cands if c.status == "ACTIVE"]
    # ── End bracketing phase ──────────────────────────────────────────────

    for it in range(1, args.iters + 1):
        if not active:
            break

        probe = max(args.probe_min, args.probe_start // (2 ** (it - 1)))

        probes_fa = tmp / f"iter{it:02d}.probes.fa"
        paf = tmp / f"iter{it:02d}.paf"

        with probes_fa.open("w") as outfa:
            for c in active:
                contig = chr2contig.get(c.chr, "")
                if not contig:
                    c.status = "NO_CONTIG"
                    continue

                if args.strategy == "bisection":
                    test_points = [(c.lo + c.hi) // 2]
                else:  # trisection
                    width = c.hi - c.lo
                    test_points = [c.lo + width // 3, c.lo + 2 * width // 3]
                c.last_test_points = test_points
                c.last_mid = test_points[0]  # legacy field

                chr_len = chr2length.get(c.chr, 0)
                for ti, t in enumerate(test_points):
                    L1 = max(0, t - probe)
                    L2 = t
                    R1 = t
                    R2 = t + probe
                    # Clamp probe extents to the chromosome so telomere-adjacent
                    # candidates keep bisecting on asymmetric probes rather than
                    # bailing out with EMPTY_PROBE.
                    if chr_len > 0:
                        L2 = min(L2, chr_len)
                        R1 = min(R1, chr_len)
                        R2 = min(R2, chr_len)

                    if (L2 - L1) < MIN_PROBE_LEN or (R2 - R1) < MIN_PROBE_LEN:
                        c.status = "EMPTY_PROBE"
                        break

                    left_seq = fasta_to_seq(faidx_region(args.child_fasta, f"{contig}:{L1}-{L2}"))
                    right_seq = fasta_to_seq(faidx_region(args.child_fasta, f"{contig}:{R1}-{R2}"))

                    if not left_seq or not right_seq:
                        c.status = "EMPTY_PROBE"
                        break

                    lid = f"{c.tag}|{c.chr}|REFINE|it{it}|t{ti}|p{t}|L|{L1}-{L2}"
                    rid = f"{c.tag}|{c.chr}|REFINE|it{it}|t{ti}|p{t}|R|{R1}-{R2}"
                    write_fasta_record(outfa, lid, left_seq)
                    write_fasta_record(outfa, rid, right_seq)

        run_minimap(args.mother_ref_fa, probes_fa, paf,
                   preset=args.preset, threads=args.threads,
                   secondary=args.secondary, cigar=args.cigar)

        scores = score_paf_by_patterns(paf, args.hap1_pattern, args.hap2_pattern)

        def _verdict(lw, rw, lc, rc, left_exp, right_exp):
            """Categorize a test point: ALL_LEFT / TRANSITION / ALL_RIGHT / SKIP."""
            if lc < args.min_conf or rc < args.min_conf:
                return "SKIP"
            if lw not in ("MATERNAL","PATERNAL") or rw not in ("MATERNAL","PATERNAL"):
                return "SKIP"
            if lw == left_exp and rw == right_exp:
                return "TRANSITION"
            if lw == left_exp and rw == left_exp:
                return "ALL_LEFT"
            if lw == right_exp and rw == right_exp:
                return "ALL_RIGHT"
            return "SKIP"

        for c in active:
            if c.status != "ACTIVE":
                continue

            test_points = c.last_test_points
            verdicts = []
            chr_len_c = chr2length.get(c.chr, 0)
            for ti, t in enumerate(test_points):
                L1 = max(0, t - probe)
                L2 = t
                R1 = t
                R2 = t + probe
                if chr_len_c > 0:
                    L2 = min(L2, chr_len_c)
                    R1 = min(R1, chr_len_c)
                    R2 = min(R2, chr_len_c)
                lid = f"{c.tag}|{c.chr}|REFINE|it{it}|t{ti}|p{t}|L|{L1}-{L2}"
                rid = f"{c.tag}|{c.chr}|REFINE|it{it}|t{ti}|p{t}|R|{R1}-{R2}"
                lw = scores.get(lid, (0,0,"MISSING",0.0))[2]
                lc = scores.get(lid, (0,0,"MISSING",0.0))[3]
                rw = scores.get(rid, (0,0,"MISSING",0.0))[2]
                rc = scores.get(rid, (0,0,"MISSING",0.0))[3]
                verdicts.append((t, lw, rw, lc, rc))

            # Use the first test point's probe outcomes for last_*/iters_done bookkeeping
            t0, lw0, rw0, lc0, rc0 = verdicts[0]
            c.iters_done = it
            c.last_left = lw0
            c.last_right = rw0
            c.last_conf = min(lc0, rc0)

            # If left_exp / right_exp are equal (NO_SWITCH inputs), auto-set from observation
            if c.left_exp == c.right_exp and c.call != "SWITCH":
                if lw0 != rw0:
                    c.left_exp, c.right_exp = lw0, rw0
                else:
                    continue

            if args.strategy == "bisection":
                t, lw, rw, lc, rc = verdicts[0]
                v = _verdict(lw, rw, lc, rc, c.left_exp, c.right_exp)
                if v == "TRANSITION":
                    c.lo = max(c.lo, t - probe // 2)
                    c.hi = min(c.hi, t + probe // 2)
                elif v == "ALL_LEFT":
                    c.lo = max(c.lo, t)
                elif v == "ALL_RIGHT":
                    c.hi = min(c.hi, t)
            else:  # trisection
                t1, lw1, rw1, lc1, rc1 = verdicts[0]
                t2, lw2, rw2, lc2, rc2 = verdicts[1]
                v1 = _verdict(lw1, rw1, lc1, rc1, c.left_exp, c.right_exp)
                v2 = _verdict(lw2, rw2, lc2, rc2, c.left_exp, c.right_exp)

                if v1 == "ALL_RIGHT":
                    # switch is to the left of t1
                    c.hi = min(c.hi, t1)
                elif v2 == "ALL_LEFT":
                    # switch is to the right of t2
                    c.lo = max(c.lo, t2)
                elif v1 == "ALL_LEFT" and v2 == "ALL_RIGHT":
                    # switch is strictly between t1 and t2 — best case, 3x reduction
                    c.lo = max(c.lo, t1)
                    c.hi = min(c.hi, t2)
                elif v1 == "ALL_LEFT" and v2 == "TRANSITION":
                    # switch is around t2
                    c.lo = max(c.lo, t2 - probe // 2)
                    c.hi = min(c.hi, t2 + probe // 2)
                elif v1 == "TRANSITION" and v2 == "ALL_RIGHT":
                    # switch is around t1
                    c.lo = max(c.lo, t1 - probe // 2)
                    c.hi = min(c.hi, t1 + probe // 2)
                elif v1 == "TRANSITION" and v2 == "TRANSITION":
                    # wide transition zone — narrow around midpoint of t1..t2
                    mid = (t1 + t2) // 2
                    c.lo = max(c.lo, mid - probe // 2)
                    c.hi = min(c.hi, mid + probe // 2)
                # else: both SKIP / inconsistent — leave interval unchanged

            if (c.hi - c.lo) <= args.resolution:
                c.status = "OK"

        active = [c for c in cands if c.status == "ACTIVE"]

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["tag","chr","bp_est_0based","bp_est_1based","bp_lo","bp_hi",
                    "iters","conf","left_exp","right_exp","last_left","last_right","status"])
        for c in cands:
            est0 = (c.lo + c.hi) // 2
            est1 = est0 + 1
            w.writerow([c.tag, c.chr, est0, est1, c.lo, c.hi,
                        c.iters_done, f"{c.last_conf:.3f}",
                        c.left_exp, c.right_exp, c.last_left, c.last_right, c.status])

    print(f"[OK] wrote {out_path}")


if __name__ == "__main__":
    main()
