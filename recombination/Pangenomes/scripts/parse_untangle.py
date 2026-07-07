#!/usr/bin/env python3
# Copyright (c) 2026, Dušan Slúka <xsluka@fi.muni.cz>
# SPDX-License-Identifier: BSD-3-Clause

"""
parse_untangle.py

Turns the BED produced by ``odgi untangle`` into a list of haplotype
crossovers in a compact switch TSV format used by the downstream
refinement step (``refine_pangenome_calls.py``).

What odgi untangle gives us
---------------------------
``odgi untangle`` walks one query path (the child's haplotype, e.g.
``PAN027#1#chr1``) along the pangenome graph and reports, for every
short segment of the query, which target paths share the most graph
nodes with it.  Each row of the BED has the columns

    query.name  query.start  query.end  target.name  target.start
    target.end  jaccard  strand  self.coverage  nth.best

The ``jaccard`` value is the Jaccard similarity coefficient
(Jaccard, 1912) between the node-set of the query segment and the
node-set of the target segment; the higher it is, the more confident
the assignment.  ``nth.best`` is 1-indexed and is used here to keep
only the single best target per segment.

Crossover-detection rule
------------------------
Segments are sorted along the chromosome and a running "current
parental haplotype" is maintained.  Whenever the next segment is
better explained by the *other* parental haplotype, a crossover is
emitted; its coordinates are the gap between the last segment of the
old haplotype (``bp_lo``) and the first segment of the new one
(``bp_hi``), with an estimated midpoint as the best single position.
The confidence is the lower of the two flanking Jaccard scores —
flagging cases where one side is well-supported but the other is not.

Query and target names follow the PanSN-spec
(https://github.com/pangenome/PanSN-spec) ``SAMPLE#HAP#CHR``
convention introduced by ``prepare_pggb_input.py``.
"""
import argparse
import collections
import re
from pathlib import Path


def parse_chr_from_pansn(name: str) -> str:
    """Extract chromosome from PanSN name like PAN010#1#chr1."""
    parts = name.split('#')
    if len(parts) >= 3:
        return parts[2]
    m = re.search(r'chr[\dXYMT]+', name, re.IGNORECASE)
    return m.group(0) if m else name


def parse_sample_hap(name: str):
    """Extract (sample, hap_number) from PanSN name like PAN010#1#chr1."""
    parts = name.split('#')
    if len(parts) >= 2:
        return parts[0], parts[1]
    return name, "?"


def main():
    ap = argparse.ArgumentParser(description="Parse odgi untangle output → crossover calls")
    ap.add_argument("--bed", required=True, help="odgi untangle BED output")
    ap.add_argument("--child-id", required=True)
    ap.add_argument("--parent-id", required=True)
    ap.add_argument("--hap", required=True, help="Child haplotype: hap1 or hap2")
    ap.add_argument("--min-identity", type=float, default=0.90,
                    help="Minimum Jaccard identity to consider a hit valid")
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    hap_num = "1" if args.hap == "hap1" else "2"
    query_prefix = f"{args.child_id}#{hap_num}#"

    # ── Read untangle BED ──────────────────────────────────────────────────
    # Keep only the rows that belong to the requested child haplotype, that
    # represent the best target match (nth.best == 1) and that clear the
    # Jaccard identity threshold.  The surviving rows are grouped by
    # chromosome, ready for a chromosome-by-chromosome sweep.
    segments = collections.defaultdict(list)

    with open(args.bed) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 10:
                continue

            q_name = parts[0]
            if not q_name.startswith(query_prefix):
                continue

            q_start = int(parts[1])
            q_end = int(parts[2])
            t_name = parts[3]
            jaccard = float(parts[6])
            nth = int(parts[9]) if parts[9].isdigit() else 1

            # Only consider best hits (nth.best == 1 in odgi untangle, 1-indexed)
            if nth != 1:
                continue

            if jaccard < args.min_identity:
                continue

            # Determine target haplotype
            t_sample, t_hap = parse_sample_hap(t_name)
            if t_sample != args.parent_id:
                continue

            q_chr = parse_chr_from_pansn(q_name)
            segments[q_chr].append((q_start, q_end, t_hap, jaccard))

    # ── Detect switches per chromosome ─────────────────────────────────────
    # Walk the segments in order; record the parental haplotype assignment
    # of each one.  A switch is reported whenever two consecutive segments
    # disagree, with the crossover bracketed by the end of the last
    # "old-haplotype" segment and the start of the first "new-haplotype"
    # segment.
    switches = []

    for chrom in sorted(segments.keys()):
        segs = sorted(segments[chrom], key=lambda x: x[0])
        if len(segs) < 2:
            continue

        prev_hap = segs[0][2]
        for i in range(1, len(segs)):
            curr_hap = segs[i][2]
            if curr_hap != prev_hap:
                bp_lo = segs[i - 1][1]   # end of last segment matching the old haplotype
                bp_hi = segs[i][0]       # start of first segment matching the new haplotype
                bp_est = (bp_lo + bp_hi) // 2  # best single-position estimate

                # Confidence: average Jaccard of flanking segments
                conf_left = segs[i - 1][3]
                conf_right = segs[i][3]
                conf = min(conf_left, conf_right)

                switches.append({
                    'tag': f"{args.child_id}_{args.hap}",
                    'chr': chrom,
                    'bp_est_0based': bp_est,
                    'bp_est_1based': bp_est + 1,
                    'bp_lo': bp_lo,
                    'bp_hi': bp_hi,
                    'conf': conf,
                    'left_hap': f"hap{prev_hap}",
                    'right_hap': f"hap{curr_hap}",
                    'method': 'pangenome',
                    'status': 'OK' if (bp_hi - bp_lo) < 100_000 else 'ACTIVE',
                })

                prev_hap = curr_hap

    # ── Write output ───────────────────────────────────────────────────────
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = "tag\tchr\tbp_est_0based\tbp_est_1based\tbp_lo\tbp_hi\tconf\tleft_hap\tright_hap\tmethod\tstatus\n"
    with out_path.open('w') as f:
        f.write(header)
        for s in switches:
            f.write('\t'.join(str(s[k]) for k in
                              ['tag', 'chr', 'bp_est_0based', 'bp_est_1based',
                               'bp_lo', 'bp_hi', 'conf', 'left_hap', 'right_hap',
                               'method', 'status']) + '\n')

    print(f"[OK] {out_path}  ({len(switches)} switches on {len(segments)} chromosomes)")


if __name__ == "__main__":
    main()
