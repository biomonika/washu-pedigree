#!/usr/bin/env python3
import argparse
import gzip
import io
import re
import sys
from typing import Dict, Tuple, List, Optional
import numpy as np
import pandas as pd
from tqdm import tqdm


# =========================
# Blocks loading & indexing
# =========================

def split_info(cell: Optional[str]):
    # "PAN010.chr1.haplotype1:142532907-142550224" -> (chrom, start, end)
    if not cell or (isinstance(cell, float) and pd.isna(cell)):
        return None, None, None
    chrom, coords = cell.split(":", 1)
    s, e = coords.split("-", 1)
    return chrom, int(s), int(e)

def load_blocks(block_paths: List[str]):
    """
    Load 1+ TSVs (no header) whose 3 columns are:
      grandparent_info, mother_info, granddaughter_info
    Returns:
      - info_by_block_and_asm[(block_id, assembly)] = (chrom, start, end)
      - chrom_index[chrom] = (starts[int64], ends[int64], block_ids[int64])  (sorted by start)
    """
    colnames = ["grandparent_info", "mother_info", "granddaughter_info"]
    frames = [pd.read_csv(p, sep="\t", header=None, names=colnames, dtype=str) for p in block_paths]
    blocks = (pd.concat(frames, ignore_index=True)
                .drop_duplicates()
                .reset_index(drop=True))

    info_by_block_and_asm: Dict[Tuple[int, str], Tuple[str, int, int]] = {}
    rows = []
    # enumerate produces a stable block_id across assemblies
    for idx, row in blocks.iterrows():
        for col in colnames:
            chrom, s, e = split_info(row[col])
            if chrom is None:
                continue
            asm = chrom.split(".", 1)[0]
            info_by_block_and_asm[(idx, asm)] = (chrom, s, e)
            rows.append((chrom, s, e, idx))

    if not rows:
        raise ValueError("No valid block intervals parsed from the provided TSVs.")

    # Build per-chrom arrays
    df = pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "BlockID"])
    # sort for binary search
    df.sort_values(["Chromosome", "Start"], inplace=True, kind="mergesort")

    chrom_index = {}
    for chrom, sub in df.groupby("Chromosome", sort=False):
        starts = sub["Start"].to_numpy(dtype=np.int64, copy=True)
        ends   = sub["End"].to_numpy(dtype=np.int64, copy=True)
        bids   = sub["BlockID"].to_numpy(dtype=np.int64, copy=True)

        # Optional safety check for overlaps within a chromosome (comment out if legit overlaps exist)
        # if np.any(starts[1:] < ends[:-1]):
        #     raise ValueError(f"Overlapping blocks detected on {chrom}; binary search assumes non-overlap.")

        chrom_index[chrom] = (starts, ends, bids)

    return info_by_block_and_asm, chrom_index


# =========================
# Fast block lookup
# =========================

def find_block_id(chrom_index, chrom: str, pos: int) -> Optional[int]:
    """
    Return BlockID for interval containing position POS, or None.
    Membership rule: Start <= pos < End  (half-open).
    """
    tup = chrom_index.get(chrom)
    if tup is None:
        return None
    starts, ends, bids = tup
    i = np.searchsorted(starts, pos, side="right") - 1
    if i >= 0 and pos < ends[i]:  # Start[i] <= pos < End[i]
        return int(bids[i])
    return None


# =========================
# Info field parsing (fast)
# =========================

def parse_info_min(info: str):
    """
    Extract QNAME, QSTART, QSTRAND from INFO without building a full dict.
    Returns (qname, qstart:int, qstrand) or (None, None, None).
    """
    qname = qstart = qstrand = None
    for item in info.split(";"):
        if not item:
            continue
        if item.startswith("QNAME="):
            qname = item[6:]
        elif item.startswith("QSTART="):
            try:
                qstart = int(item[7:])
            except ValueError:
                qstart = None
        elif item.startswith("QSTRAND="):
            qstrand = item[8:]
    return qname, qstart, qstrand


def parse_qname_as_block(qname: str) -> Tuple[str, int, int]:
    # Accept "PAN028.chr1.haplotype2:123-456" or "chr1:123-456"
    m = re.match(r"([^:]+):(\d+)-(\d+)$", qname)
    if not m:
        raise ValueError(f"Invalid QNAME format: {qname}. Supply --blocks TSVs explicitly or ensure 'chr:start-end'.")
    chrom, s, e = m.group(1), int(m.group(2)), int(m.group(3))
    return chrom, s, e


# =========================
# I/O helpers
# =========================

def open_text_maybe_gzip(path: str, mode: str):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, mode.replace("t", "")), encoding="utf-8")
    return open(path, mode, encoding="utf-8", newline="")

def count_non_header_lines(path: str) -> Optional[int]:
    if path in ("-", None):
        return None
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    with opener(path, "rt", encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                n += 1
    return n


# =========================
# Core conversion
# =========================

def convert_vcf(input_vcf: str,
                output_vcf: str,
                block_paths: Optional[List[str]] = None,
                skip_missing: bool = False,
                show_progress: bool = True,
                count_lines_first: bool = False):

    if block_paths:
        info_by_blk_asm, chrom_index = load_blocks(block_paths)
    else:
        info_by_blk_asm, chrom_index = {}, {}

    total = count_non_header_lines(input_vcf) if (show_progress and count_lines_first) else None
    pbar = tqdm(total=total, unit="var", desc="Converting VCF", disable=not show_progress, mininterval=0.2)

    fin = open_text_maybe_gzip(input_vcf, "rt")
    fout = open_text_maybe_gzip(output_vcf, "wt")
    try:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                if not skip_missing:
                    fout.write(line)
                if show_progress: pbar.update(1)
                continue

            chrom, pos_s = fields[0], fields[1]
            pos1 = int(pos_s)
            vid, ref, alt, qual, flt, info = fields[2:8]
            fmt = fields[8] if len(fields) > 8 else ""
            sample_and_more = fields[9:] if len(fields) > 9 else []

            qname, qstart, qstrand = parse_info_min(info)
            if qname is None or qstart is None or qstrand is None:
                if not skip_missing:
                    fout.write(line)
                if show_progress: pbar.update(1)
                continue

            # Get block coordinates for destination assembly
            if chrom_index:
                block_id = find_block_id(chrom_index, chrom, pos1)
                if block_id is None:
                    if skip_missing:
                        if show_progress: pbar.update(1)
                        continue
                    else:
                        fout.write(line)
                        if show_progress: pbar.update(1)
                        continue
                dest_asm = qname.split(".", 1)[0]
                dest = info_by_blk_asm.get((block_id, dest_asm))
                if dest is None:
                    if skip_missing:
                        if show_progress: pbar.update(1)
                        continue
                    else:
                        fout.write(line)
                        if show_progress: pbar.update(1)
                        continue
                new_chrom, aln_start, aln_end = dest
            else:
                # fallback: parse coordinates right from QNAME
                new_chrom, aln_start, aln_end = parse_qname_as_block(qname)

            # Map position based on strand
            if qstrand == "+":
                new_pos = aln_start + qstart
            elif qstrand == "-":
                new_pos = aln_end - qstart
            else:
                # unknown strand; keep record or skip
                if not skip_missing:
                    fout.write(line)
                if show_progress: pbar.update(1)
                continue

            # Rebuild INFO (drop QNAME/QSTART/QSTRAND, add ORIG)
            # Fast rebuild: scan again, skipping those keys
            orig_pos = f"{chrom}:{pos_s}"
            kept = []
            for item in info.split(";"):
                if not item or item.startswith(("QNAME=", "QSTART=", "QSTRAND=")):
                    continue
                kept.append(item)
            kept.append(f"ORIG={orig_pos}")
            new_info = ";".join(kept)

            # Switch REF and ALT
            ref, alt = alt, ref

            out = [new_chrom, str(new_pos), vid, ref, alt, qual, flt, new_info]
            if fmt:
                out.append(fmt)
            out.extend(sample_and_more)
            fout.write("\t".join(out) + "\n")

            if show_progress: pbar.update(1)

    finally:
        if fin is not sys.stdin:
            fin.close()
        if fout is not sys.stdout:
            fout.close()
        if show_progress:
            pbar.close()


# =========================
# CLI
# =========================

def main():
    ap = argparse.ArgumentParser(
        description="Fast map of VCF query coordinates to reference using shared block tables."
    )
    ap.add_argument("input_vcf", help="Input VCF/VCF.GZ (use '-' for stdin)")
    ap.add_argument("output_vcf", help="Output VCF (use '-' for stdout)")
    ap.add_argument("-b", "--blocks", action="append", default=[],
                    help="Repeatable: -b blocks.tsv (3 columns: grandparent/mother/granddaughter).")
    ap.add_argument("--skip-missing", action="store_true",
                    help="Skip variants that don't map to a block (default: keep unchanged).")
    ap.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bar.")
    ap.add_argument("--count-lines-first", action="store_true",
                    help="Count non-header lines first for determinate progress (slower for .gz).")
    args = ap.parse_args()

    convert_vcf(
        input_vcf=args.input_vcf,
        output_vcf=args.output_vcf,
        block_paths=args.blocks or None,
        skip_missing=args.skip_missing,
        show_progress=not args.no_progress,
        count_lines_first=args.count_lines_first
    )

if __name__ == "__main__":
    main()