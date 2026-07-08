#!/usr/bin/env python3
# Copyright (c) 2026, Dušan Slúka <xsluka@fi.muni.cz>
# SPDX-License-Identifier: BSD-3-Clause

"""Karyogram visualisation of pipeline crossover/breakpoint runs.

Renders one PNG per run, showing all autosomes (chr1..chr22) with the
maternal and paternal ideograms placed side-by-side in a mirrored layout:

        [-- maternal --]   chrN   [-- paternal --]
         (p-tel outer)            (p-tel outer)

Each chromosome bar is coloured by inter-breakpoint segments alternating
between hap1 (green) and hap2 (blue). Below the bar, the centromere
is drawn as a thin dark underscore and small Mb tick labels indicate
approximate position along the chromosome. Red vertical lines mark
breakpoint positions (alpha proportional to confidence; a pale red box
spans the BED uncertainty interval).

Expected input:
    A pipeline run directory containing crossovers.bed
    (chrom start end name score strand, where `name` is one of
    mother_hap{1,2}_SWITCH / father_hap{1,2}_SWITCH and `score` = conf*1000).
    The script auto-detects <run_dir>/crossovers.bed and
    <run_dir>/pangenome/crossovers.bed.

Cytoband file (optional):
    UCSC hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
    Used only to extract centromere position (acen stain). The hap-coloured
    ideogram replaces G-banding, so cytoband bands themselves are not drawn.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


CHROM_ORDER = [f"chr{i}" for i in range(1, 23)]

CHR_NAME_RE = re.compile(r"\.(chr0*[0-9XY]+)\.")
CHR_NORM_RE = re.compile(r"^chr0*([1-9]|1[0-9]|2[0-2]|X|Y)$")


def norm_chrom(chrom: str) -> str | None:
    """Canonical 'chrN' (no zero-padding) for any of chr1/chr01/chrX/etc."""
    m = CHR_NORM_RE.match(chrom.strip())
    return f"chr{m.group(1)}" if m else None


C_HAP1       = "#2ca02c"
C_HAP2       = "#1f77b4"
C_BREAKPOINT = "#e62828"
C_CENTROMERE = "#2a2a2a"
C_BAR_EDGE   = "#222222"
C_TICK       = "#888888"
C_TICK_LABEL = "#555555"

HG38_CENTROMERES = {
    "chr1":  (122026459, 125184587),
    "chr2":  (92188145,  94090557),
    "chr3":  (90772458,  93655574),
    "chr4":  (49708101,  51743951),
    "chr5":  (46485900,  50059807),
    "chr6":  (58553888,  59829934),
    "chr7":  (58169653,  60828234),
    "chr8":  (44033744,  45877265),
    "chr9":  (43236167,  45518558),
    "chr10": (39686682,  41593521),
    "chr11": (51078348,  54425074),
    "chr12": (34769407,  37185252),
    "chr13": (16000000,  18051248),
    "chr14": (16000000,  18173523),
    "chr15": (17083673,  19725254),
    "chr16": (36311158,  38265669),
    "chr17": (22813679,  26616164),
    "chr18": (15460898,  20861206),
    "chr19": (24631782,  27190874),
    "chr20": (26369569,  28464347),
    "chr21": (10864560,  12915808),
    "chr22": (13684672,  15054318),
}


def parse_fai(fai: Path) -> dict[str, int]:
    """Reads a samtools-style FASTA index and returns ``{chrN: length}``.

    Contig names may be assembly-decorated (``PAN027.chr1.maternal``),
    PanSN-formatted (``PAN027#1#chr1``), or plain (``chr1``); all three
    forms are normalised to ``chrN``.  When several haplotypes report
    the same chromosome the longest length is kept so the karyogram
    bars never get truncated."""
    out: dict[str, int] = {}
    for line in fai.read_text().splitlines():
        if not line.strip():
            continue
        name, length = line.split("\t")[:2]
        raw = None
        m = CHR_NAME_RE.search(name)
        if m:
            raw = m.group(1)
        elif "#" in name:
            raw = name.rsplit("#", 1)[-1]
        elif name.startswith("chr"):
            raw = name
        chrom = norm_chrom(raw) if raw else None
        if chrom and chrom in CHROM_ORDER:
            out[chrom] = max(out.get(chrom, 0), int(length))
    return out


def load_centromeres_from_cytoband(path: Path | None) -> dict[str, tuple[int, int]]:
    """Extract centromere region (acen stain) per chromosome from UCSC cytoBand."""
    out: dict[str, tuple[int, int]] = {}
    if path is None:
        return out
    opener = gzip.open if str(path).endswith(".gz") else open
    acens: dict[str, list[tuple[int, int]]] = {}
    with opener(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or p[4] != "acen":
                continue
            chrom = norm_chrom(p[0])
            if chrom not in CHROM_ORDER:
                continue
            acens.setdefault(chrom, []).append((int(p[1]), int(p[2])))
    for c, regs in acens.items():
        out[c] = (min(s for s, _ in regs), max(e for _, e in regs))
    return out


def _extract_chrom_and_hap(raw: str) -> tuple[str | None, str | None]:
    """Parse a BED chr column. Returns (canonical chrom, haplotype) where
    haplotype is 'maternal'/'paternal' if the name is decorated (e.g.
    'PAN027.chr1.maternal'), else None."""
    raw = raw.strip()
    m = CHR_NAME_RE.search(raw)
    if m:
        chrom = norm_chrom(m.group(1))
        hap = None
        if raw.endswith(".maternal"):
            hap = "maternal"
        elif raw.endswith(".paternal"):
            hap = "paternal"
        return chrom, hap
    return norm_chrom(raw), None


def load_centromeres_bed(path: Path | None, stain_prefix: str = "active_hor"
                         ) -> dict[str, dict[str, tuple[int, int]]]:
    """Load centromere intervals from a cenSat-style BED. Returns
    {chrom: {'maternal': (s,e), 'paternal': (s,e), 'both': (s,e)}}.

    - If the BED chr column is decorated with haplotype (e.g.
      'PAN027.chrN.maternal'), per-haplotype spans are kept separately
      so the maternal and paternal panels show their own centromere.
    - If the BED is plain (e.g. 'chr1'), the same span is reused for
      both haplotypes under the 'both' key.
    - Only rows whose column-4 label starts with `stain_prefix` are kept
      (default 'active_hor' = the modern T2T cenSat centromere definition).
      Pass stain_prefix='' to keep all rows (legacy/plain BEDs).
    - Multiple intervals per (chrom, hap) are merged to (min start, max end).
    """
    out: dict[str, dict[str, tuple[int, int]]] = {}
    if path is None:
        return out
    opener = gzip.open if str(path).endswith(".gz") else open
    spans: dict[tuple[str, str], list[tuple[int, int]]] = {}
    with opener(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            p = line.split("\t")
            if len(p) < 3:
                continue
            chrom, hap = _extract_chrom_and_hap(p[0])
            if chrom not in CHROM_ORDER:
                continue
            if stain_prefix and (len(p) < 4 or not p[3].startswith(stain_prefix)):
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            spans.setdefault((chrom, hap or "both"), []).append((s, e))
    for (chrom, hap), regs in spans.items():
        out.setdefault(chrom, {})[hap] = (
            min(s for s, _ in regs), max(e for _, e in regs))
    return out


def load_breakpoints(bed: Path) -> dict[str, dict[str, list[tuple[int, int, float]]]]:
    """Reads a pipeline ``crossovers.bed`` and groups the SWITCH calls by
    parental branch (maternal / paternal) and chromosome.

    The BED ``name`` column must follow ``mother_hapN_SWITCH`` /
    ``father_hapN_SWITCH``; ``mother_*`` rows feed the maternal panel
    and ``father_*`` rows feed the paternal panel.  The ``score``
    column carries the confidence multiplied by 1000 (the controller
    encodes confidence as an integer to stay BED-compatible) and is
    decoded back into the ``[0, 1]`` range.  Rows without ``SWITCH``
    in the name are skipped so NO_SWITCH / AMBIG / CROSS_MAP calls do
    not show up in the figure."""
    out: dict[str, dict[str, list]] = {"maternal": {}, "paternal": {}}
    with bed.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 4:
                continue
            chrom = norm_chrom(p[0])
            if chrom not in CHROM_ORDER:
                continue
            try:
                start, end = int(p[1]), int(p[2])
            except ValueError:
                continue
            name = p[3]
            if "SWITCH" not in name:
                continue
            try:
                score = int(p[4]) if len(p) >= 5 else 0
            except ValueError:
                score = 0
            conf = max(0.0, min(1.0, score / 1000.0))
            if name.startswith("mother_"):
                ch = "maternal"
            elif name.startswith("father_"):
                ch = "paternal"
            else:
                continue
            out[ch].setdefault(chrom, []).append((start, end, conf))
    for ch in out.values():
        for c in ch:
            ch[c].sort()
    return out


def x_mat(X_bp: int, L_bp: int, gap_mb: float) -> float:
    """Maternal panel: p-tel outer-left, q-tel inner. Reads naturally L→R."""
    return -gap_mb - (L_bp - X_bp) / 1e6


def x_pat(X_bp: int, L_bp: int, gap_mb: float) -> float:
    """Paternal panel: mirror of maternal. p-tel outer-right, q-tel inner."""
    return +gap_mb + (L_bp - X_bp) / 1e6


def pick_tick_step_mb(L_mb: float) -> int:
    """Chooses a Mb tick spacing that keeps roughly 3–6 labels per
    chromosome regardless of its length."""
    if L_mb <= 60:
        return 20
    if L_mb <= 130:
        return 25
    return 50


def draw_panel(ax, side: str, y: float, L_bp: int,
               cent: tuple[int, int] | None,
               breakpoints: list[tuple[int, int, float]],
               height: float, gap_mb: float) -> None:
    """Draws one side (maternal or paternal) of a chromosome row.

    The panel is composed of three visual layers:
      * the chromosome bar itself, painted in alternating hap1/hap2
        colours between consecutive breakpoint midpoints,
      * a thin axis line beneath it, with a dark centromere block
        sitting on the line and Mb tick labels under it,
      * the breakpoint markers — a translucent uncertainty rectangle
        spanning ``bp_start..bp_end`` plus a sharp midpoint line, both
        with opacity proportional to the call's confidence.

    The horizontal projection is mirrored between the two sides so the
    maternal panel grows leftward and the paternal panel grows
    rightward from the central chromosome label."""
    proj = x_mat if side == "maternal" else x_pat
    L_mb = L_bp / 1e6
    x_inner = proj(L_bp, L_bp, gap_mb)
    x_outer = proj(0, L_bp, gap_mb)
    x_left = min(x_inner, x_outer)

    # Alternating hap1/hap2 segments INSIDE the chromosome bar
    bp_mids = sorted({(s + e) // 2 for (s, e, _) in breakpoints})
    edges = [0] + bp_mids + [L_bp]
    for i in range(len(edges) - 1):
        s, e = edges[i], edges[i + 1]
        if e <= s:
            continue
        x1, x2 = proj(s, L_bp, gap_mb), proj(e, L_bp, gap_mb)
        xb, w = min(x1, x2), abs(x2 - x1)
        col = C_HAP1 if i % 2 == 0 else C_HAP2
        ax.add_patch(mpatches.Rectangle(
            (xb, y - height / 2), w, height,
            facecolor=col, edgecolor="none", zorder=2))

    # Chromosome outline
    ax.add_patch(mpatches.Rectangle(
        (x_left, y - height / 2), L_mb, height,
        facecolor="none", edgecolor=C_BAR_EDGE, linewidth=0.7, zorder=4))

    # Horizontal gray axis line below the bar — spans the chr length
    axis_y = y - height / 2 - 0.16
    ax.plot([x_left, x_left + L_mb], [axis_y, axis_y],
            color=C_TICK, linewidth=0.8, zorder=2,
            solid_capstyle="butt")

    # Centromere as a dark rectangle SITTING ON the axis line
    if cent:
        cs, ce = cent
        x1, x2 = proj(cs, L_bp, gap_mb), proj(ce, L_bp, gap_mb)
        xb, w_cent = min(x1, x2), abs(x2 - x1)
        # Minimum visible width: 2 Mb (so it never disappears on long chroms)
        if w_cent < 2.0:
            x_center = (x1 + x2) / 2
            xb = x_center - 1.0
            w_cent = 2.0
        ax.add_patch(mpatches.Rectangle(
            (xb, axis_y - 0.04), w_cent, 0.08,
            facecolor=C_CENTROMERE, edgecolor="none", zorder=3))

    # Mb tick marks (drop below the axis line) + bold labels under them
    step = pick_tick_step_mb(L_mb)
    last_drawn = -1
    pos = 0
    while pos <= L_mb + 1e-6:
        x_t = proj(int(pos * 1e6), L_bp, gap_mb)
        ax.plot([x_t, x_t], [axis_y - 0.04, axis_y - 0.07],
                color=C_TICK, linewidth=0.7, zorder=2,
                solid_capstyle="butt")
        ax.text(x_t, axis_y - 0.09, f"{int(pos)}",
                ha="center", va="top", fontsize=10, fontweight="bold",
                color=C_TICK_LABEL, zorder=2)
        last_drawn = pos
        pos += step
    # Show the q-telomere length too if it's not too close to the last tick
    last_mb = int(L_mb)
    if last_mb - last_drawn > step * 0.4:
        x_t = proj(L_bp, L_bp, gap_mb)
        ax.plot([x_t, x_t], [axis_y - 0.04, axis_y - 0.07],
                color=C_TICK, linewidth=0.7, zorder=2,
                solid_capstyle="butt")
        ax.text(x_t, axis_y - 0.09, f"{last_mb}",
                ha="center", va="top", fontsize=10, fontweight="bold",
                color=C_TICK_LABEL, zorder=2)

    # Red breakpoint markers (uncertainty rect + sharp midpoint line)
    bar_top = y + height / 2
    bar_bot = y - height / 2
    for (s, e, conf) in breakpoints:
        x1, x2 = proj(s, L_bp, gap_mb), proj(e, L_bp, gap_mb)
        xb, w = min(x1, x2), abs(x2 - x1)
        if w > 0:
            ax.add_patch(mpatches.Rectangle(
                (xb, bar_bot - 0.02), w, (bar_top - bar_bot) + 0.04,
                facecolor=C_BREAKPOINT, edgecolor="none",
                alpha=0.18 + 0.22 * conf, zorder=5))
        x_bp = proj((s + e) // 2, L_bp, gap_mb)
        ax.plot([x_bp, x_bp],
                [bar_bot - 0.04, bar_top + 0.04],
                color=C_BREAKPOINT, linewidth=1.6,
                alpha=0.65 + 0.35 * conf, zorder=6,
                solid_capstyle="butt")


def _pick_cent(cent_from_bed: dict, cyto_cents: dict, chrom: str,
               hap: str) -> tuple[int, int] | None:
    """Resolution order: cenSat BED for the matching haplotype → cenSat 'both'
    (if BED was plain) → cytoBand acen → hardcoded hg38 fallback."""
    per_chr = cent_from_bed.get(chrom)
    if per_chr:
        if hap in per_chr:
            return per_chr[hap]
        if "both" in per_chr:
            return per_chr["both"]
        # cenSat had only the other hap — fall back to whatever single span
        # is available (better than nothing)
        return next(iter(per_chr.values()))
    return cyto_cents.get(chrom) or HG38_CENTROMERES.get(chrom)


def render(run_dir: Path, bed: Path, fai: Path, cyto: Path | None,
           cent_bed: Path | None, out: Path, title: str | None) -> None:
    """Renders a karyogram of one pipeline run and writes the PNG to
    ``out``.

    The figure stacks the 22 autosomes top-to-bottom, each as a single
    row with the maternal panel on the left, the chromosome label in
    the centre, and the paternal panel on the right.  Per-row layout
    parameters (``row_h``, ``bar_h``, ``gap_mb``) are tuned so the
    output stays legible at the default DPI for chromosomes up to
    ~250 Mb."""
    lengths = parse_fai(fai)
    cent_from_bed = load_centromeres_bed(cent_bed)
    cyto_cents = load_centromeres_from_cytoband(cyto)
    calls = load_breakpoints(bed)

    chroms = [c for c in CHROM_ORDER if c in lengths]
    if not chroms:
        raise SystemExit(f"No autosome lengths parsed from {fai}")

    max_len_mb = max(lengths[c] for c in chroms) / 1e6
    gap_mb = 6.0
    n = len(chroms)
    row_h = 1.30
    bar_h = 0.34

    fig_w = min(22.0, max(14.0, 0.05 * (2 * max_len_mb + 2 * gap_mb)))
    fig_h = row_h * n + 2.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for i, chrom in enumerate(chroms):
        y = (n - 1 - i) * row_h
        L_bp = lengths[chrom]
        cent_mat = _pick_cent(cent_from_bed, cyto_cents, chrom, "maternal")
        cent_pat = _pick_cent(cent_from_bed, cyto_cents, chrom, "paternal")
        draw_panel(ax, "maternal", y, L_bp, cent_mat,
                   calls["maternal"].get(chrom, []), bar_h, gap_mb)
        draw_panel(ax, "paternal", y, L_bp, cent_pat,
                   calls["paternal"].get(chrom, []), bar_h, gap_mb)
        n_digits = "".join(ch for ch in chrom if ch.isdigit())
        label = f"chr{int(n_digits):02d}" if n_digits else chrom
        ax.text(0.0, y, label, ha="center", va="center",
                fontsize=12, fontweight="bold", zorder=10)

    x_pad = 5.0
    ax.set_xlim(-(max_len_mb + gap_mb + x_pad),
                +(max_len_mb + gap_mb + x_pad))
    ax.set_ylim(-row_h, n * row_h)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("maternal           chr           paternal", fontsize=12)

    n_mat = sum(len(v) for v in calls["maternal"].values())
    n_pat = sum(len(v) for v in calls["paternal"].values())
    ax.set_title(
        title or (f"{run_dir.name}    "
                  f"maternal={n_mat}  paternal={n_pat}  total={n_mat + n_pat}"),
        fontsize=14, fontweight="bold")

    handles = [
        mpatches.Patch(facecolor=C_HAP1,
                       label="hap1 segment (between breakpoints)"),
        mpatches.Patch(facecolor=C_HAP2,
                       label="hap2 segment (between breakpoints)"),
        plt.Line2D([0], [0], color=C_BREAKPOINT, lw=2.0,
                   label="breakpoint  (α ∝ confidence)"),
        mpatches.Patch(facecolor=C_CENTROMERE,
                       label=("centromere (cenSat active_hor)"
                              if cent_from_bed else
                              "centromere (hg38 UCSC chromInfo)")),
    ]
    ax.legend(handles=handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.005),
              ncol=4, fontsize=11, frameon=True,
              edgecolor="#bbbbbb", framealpha=0.95,
              handlelength=2.5, handletextpad=0.8,
              columnspacing=1.8)

    out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[wrote] {out}  ({n_mat + n_pat} breakpoints)")


def autodetect_bed(run_dir: Path) -> Path | None:
    """Finds the ``crossovers.bed`` produced by a pipeline run.

    Two conventional locations are tried:
      1. ``run_dir/crossovers.bed`` (root of the run directory)
      2. ``run_dir/pangenome/crossovers.bed`` (pangenome sub-run)

    Both are checked so the plotter works regardless of the exact
    directory layout the caller used."""
    for cand in (run_dir / "crossovers.bed",
                 run_dir / "pangenome" / "crossovers.bed"):
        if cand.exists():
            return cand
    return None


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, action="append",
                    help="Pipeline run directory. Repeat for multiple runs.")
    ap.add_argument("--fai", required=True,
                    help=".fai with chromosome lengths "
                         "(decorated names like PAN027.chrN.maternal are ok).")
    ap.add_argument("--cytobands",
                    help="UCSC cytoBand.txt[.gz]. Only used to extract "
                         "centromere position; otherwise hardcoded hg38 "
                         "centromeres are used.")
    ap.add_argument("--centromeres",
                    help="BED of centromere intervals "
                         "(e.g. T2T-CHM13 cenSat active_hor track). "
                         "Takes precedence over --cytobands and hg38 fallback.")
    ap.add_argument("--bed",
                    help="Override BED path (single --run-dir only).")
    ap.add_argument("--out",
                    help="Output PNG path (single --run-dir only).")
    ap.add_argument("--out-dir",
                    help="Output dir for multi-run mode "
                         "(default: <run_dir>/figures/).")
    ap.add_argument("--title", help="Title override.")
    args = ap.parse_args()

    fai = Path(args.fai)
    if not fai.exists():
        raise SystemExit(f"FAI not found: {fai}")
    cyto = Path(args.cytobands) if args.cytobands else None
    if cyto and not cyto.exists():
        raise SystemExit(f"Cytoband file not found: {cyto}")
    cent_bed = Path(args.centromeres) if args.centromeres else None
    if cent_bed and not cent_bed.exists():
        raise SystemExit(f"Centromere BED not found: {cent_bed}")

    for rd in args.run_dir:
        run_dir = Path(rd).resolve()
        if not run_dir.is_dir():
            raise SystemExit(f"Run dir not found: {run_dir}")
        if args.bed and len(args.run_dir) == 1:
            bed = Path(args.bed)
        else:
            bed = autodetect_bed(run_dir)
            if bed is None:
                raise SystemExit(
                    f"crossovers.bed not found in {run_dir} "
                    f"(checked: ./crossovers.bed, ./pangenome/crossovers.bed)")
        if args.out and len(args.run_dir) == 1:
            out = Path(args.out)
        else:
            out_dir = Path(args.out_dir) if args.out_dir else (run_dir / "figures")
            out = out_dir / f"karyogram_{run_dir.name}.png"

        render(run_dir, bed, fai, cyto, cent_bed, out, args.title)


if __name__ == "__main__":
    main()
