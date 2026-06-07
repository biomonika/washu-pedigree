import sys
import os
import bisect
from collections import defaultdict

# ============================================================
# USER CONFIGURATION
# ============================================================
BLOCKS_DIR_MAT  = "/path/to/maternal_transmitted_blocks_per_chr"
BLOCKS_DIR_PAT  = "/path/to/paternal_transmitted_blocks_per_chr"
ANNOTATIONS_DIR = "/path/to/annotations"
PROBLEMATIC_BED = f"{ANNOTATIONS_DIR}/problematic.PAN027.bed"
SW_PRONE_BED    = f"{ANNOTATIONS_DIR}/sw_prone_regions.bed"
TRF_BED         = f"{ANNOTATIONS_DIR}/PAN027_HiFi_element_polished.v1.0.unit1_6.trf.bed"
VCF_PATTERN     = "PAN027.chr{chr}.project_{ref}.wave.filtered.vcf"  # adjust if your files are named differently
# ============================================================

# --- ARGUMENTS ---
# Usage: python3 make_stats.per_chr.py <chr> <hap> [out_base]
#   chr      : chromosome number or X (e.g. "1", "X")
#   hap      : "mat" or "pat"
#   out_base : results folder (default: results_filtered_all)
if len(sys.argv) >= 3:
    _chr, _hap = sys.argv[1], sys.argv[2]
else:
    _chr, _hap = "1", "pat"
_out_base = sys.argv[3] if len(sys.argv) >= 4 else "results_filtered_all"

_ref = "PAN027_mat" if _hap == "mat" else "PAN027_pat"
_vcf_name = VCF_PATTERN.format(chr=_chr, ref=_ref)

vcf_file        = f"{_out_base}/chr{_chr}_{_hap}/{_vcf_name}"
output_stats_tsv = f"{_out_base}/chr{_chr}_{_hap}/mutation_statistics.tsv"
output_lengths_tsv = f"{_out_base}/chr{_chr}_{_hap}/indel_lengths.tsv"
output_bed      = f"{_out_base}/chr{_chr}_{_hap}/heatmap_positions.bed"

os.makedirs(os.path.dirname(output_stats_tsv), exist_ok=True)
os.makedirs(os.path.dirname(output_bed), exist_ok=True)

# --- BLOCK LENGTHS ---
blocks_dir  = BLOCKS_DIR_MAT if _hap == "mat" else BLOCKS_DIR_PAT
blocks_file = f"{blocks_dir}/chr{_chr}.tsv"

block_intervals = []
total_block_length = 0
with open(blocks_file) as bf:
    for line in bf:
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        chrom_b, reg = parts[1].split(':')
        s, e = reg.split('-')
        s, e = int(s), int(e)
        block_intervals.append((chrom_b, s, e))
        total_block_length += e - s

# Exclusion regions for effective block length calculation
exclude_by_chrom = defaultdict(list)
for bed_path in [PROBLEMATIC_BED, SW_PRONE_BED]:
    try:
        with open(bed_path) as ef:
            for line in ef:
                if line.startswith('#') or not line.strip():
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 3:
                    continue
                exclude_by_chrom[cols[0]].append((int(cols[1]), int(cols[2])))
    except FileNotFoundError:
        pass

def subtract_intervals(block_ivs, exclude_by_chrom):
    total = 0
    for (chrom, bs, be) in block_ivs:
        remaining = [(bs, be)]
        for (es, ee) in exclude_by_chrom.get(chrom, []):
            new_remaining = []
            for (rs, re) in remaining:
                if ee <= rs or es >= re:
                    new_remaining.append((rs, re))
                else:
                    if rs < es:
                        new_remaining.append((rs, es))
                    if ee < re:
                        new_remaining.append((ee, re))
            remaining = new_remaining
        total += sum(e - s for s, e in remaining)
    return total

effective_block_length = subtract_intervals(block_intervals, exclude_by_chrom)

# --- LOAD TRF BED (tandem repeat annotations, PAN027 coordinates) ---
_trf_raw = defaultdict(list)
try:
    with open(TRF_BED) as tf:
        for line in tf:
            cols = line.strip().split('\t')
            if len(cols) < 3:
                continue
            _trf_raw[cols[0]].append((int(cols[1]), int(cols[2])))
except FileNotFoundError:
    pass

# Pre-process TRF intervals for O(log n) overlap lookup
_trf = {}
for contig, ivs in _trf_raw.items():
    ivs.sort()
    starts = [s for s, e in ivs]
    max_ends = []
    cur_max = 0
    for _, e in ivs:
        cur_max = max(cur_max, e)
        max_ends.append(cur_max)
    _trf[contig] = (starts, max_ends)

def in_trf(chrom, var_start, var_end):
    """Return True if [var_start, var_end) overlaps any TRF interval."""
    if chrom not in _trf:
        return False
    starts, max_ends = _trf[chrom]
    i = bisect.bisect_left(starts, var_end) - 1
    if i < 0:
        return False
    return max_ends[i] > var_start

# --- PARSE VCF AND CLASSIFY VARIANTS ---
snv_count = 0
indel_1bp_ins_count = 0
indel_1bp_del_count = 0
indel_2_49bp_ins_count = 0
indel_2_49bp_del_count = 0
sv_ins_count = 0
sv_del_count = 0
inv_count = 0
trf_count = 0
indel_lengths = []

with open(vcf_file, 'r') as f, open(output_bed, 'w') as bed_f:
    for line in f:
        if line.startswith('#'):
            continue

        columns = line.strip().split('\t')
        if len(columns) < 5:
            continue

        chrom = columns[0]
        pos   = int(columns[1])
        ref   = columns[3]
        alt   = columns[4].split(',')[0]
        info  = columns[7] if len(columns) > 7 else ""

        # Parse TYPE and LEN from INFO field (set by vcfwave)
        info_type = None
        info_len  = None
        for field in info.split(';'):
            if field.startswith('TYPE='):
                v = field[5:].split(',')[0].lower()  # take first value (multiallelic safety)
                if v and v != '.':
                    info_type = v
            elif field.startswith('LEN='):
                try:
                    info_len = int(field[4:].split(',')[0])  # take first value (multiallelic safety)
                except ValueError:
                    pass

        # Fallback: derive LEN from sequences when LEN is missing or 0
        len_diff = abs(len(ref) - len(alt))
        if info_len is None or info_len == 0:
            info_len = len_diff

        # Fallback: derive TYPE from sequences when TYPE is missing
        if info_type is None:
            if alt.upper().startswith('<INV') or 'SVTYPE=INV' in info:
                info_type = 'inv'
            elif len(ref) == 1 and len(alt) == 1:
                info_type = 'snp'
            elif len(alt) > len(ref):
                info_type = 'ins'
            else:
                info_type = 'del'

        mut_type = "UNKNOWN"

        if info_type == 'inv':
            inv_count += 1
            mut_type = "INV"
        elif info_type == 'snp':
            snv_count += 1
            mut_type = "SNV"
        elif info_type == 'mnp':
            snv_count += info_len  # each base in an MNP counts as one SNV
            mut_type = "SNV"
        elif info_type in ('ins', 'del'):
            if info_len == 1:
                indel_lengths.append(info_len)
                if info_type == 'ins':
                    indel_1bp_ins_count += 1
                    mut_type = "INS_1bp"
                else:
                    indel_1bp_del_count += 1
                    mut_type = "DEL_1bp"
            elif 2 <= info_len <= 49:
                indel_lengths.append(info_len)
                if info_type == 'ins':
                    indel_2_49bp_ins_count += 1
                    mut_type = "INS"
                else:
                    indel_2_49bp_del_count += 1
                    mut_type = "DEL"
            else:  # >= 50 bp = SV
                indel_lengths.append(info_len)
                if info_type == 'ins':
                    sv_ins_count += 1
                    mut_type = "SV_INS"
                else:
                    sv_del_count += 1
                    mut_type = "SV_DEL"

        # Convert VCF 1-based to BED 0-based half-open
        start = pos - 1
        end   = start + len(ref)

        if mut_type != "UNKNOWN":
            bed_f.write(f"{chrom}\t{start}\t{end}\t{mut_type}\n")
            if in_trf(chrom, start, end):
                trf_count += 1
                bed_f.write(f"{chrom}\t{start}\t{end}\tTRF\n")

# --- WRITE STATISTICS TSV ---
indel_1bp_count    = indel_1bp_ins_count + indel_1bp_del_count
indel_2_49bp_count = indel_2_49bp_ins_count + indel_2_49bp_del_count
sv_count           = sv_ins_count + sv_del_count
total_variants     = snv_count + indel_1bp_count + indel_2_49bp_count + sv_count

def rate(count, length):
    return count / length if length > 0 else 0

with open(output_stats_tsv, 'w') as out_f:
    out_f.write("Mutation_Type\tCount\tRate_per_bp_total\tRate_per_bp_effective\n")
    out_f.write(f"SNV\t{snv_count}\t{rate(snv_count, total_block_length):.6e}\t{rate(snv_count, effective_block_length):.6e}\n")
    out_f.write(f"Indel_1bp_Total\t{indel_1bp_count}\t{rate(indel_1bp_count, total_block_length):.6e}\t{rate(indel_1bp_count, effective_block_length):.6e}\n")
    out_f.write(f"Indel_1bp_Insertion\t{indel_1bp_ins_count}\t{rate(indel_1bp_ins_count, total_block_length):.6e}\t{rate(indel_1bp_ins_count, effective_block_length):.6e}\n")
    out_f.write(f"Indel_1bp_Deletion\t{indel_1bp_del_count}\t{rate(indel_1bp_del_count, total_block_length):.6e}\t{rate(indel_1bp_del_count, effective_block_length):.6e}\n")
    out_f.write(f"Indel_2_49bp_Total\t{indel_2_49bp_count}\t{rate(indel_2_49bp_count, total_block_length):.6e}\t{rate(indel_2_49bp_count, effective_block_length):.6e}\n")
    out_f.write(f"Indel_2_49bp_Insertion\t{indel_2_49bp_ins_count}\t{rate(indel_2_49bp_ins_count, total_block_length):.6e}\t{rate(indel_2_49bp_ins_count, effective_block_length):.6e}\n")
    out_f.write(f"Indel_2_49bp_Deletion\t{indel_2_49bp_del_count}\t{rate(indel_2_49bp_del_count, total_block_length):.6e}\t{rate(indel_2_49bp_del_count, effective_block_length):.6e}\n")
    out_f.write(f"SV_Total\t{sv_count}\t{rate(sv_count, total_block_length):.6e}\t{rate(sv_count, effective_block_length):.6e}\n")
    out_f.write(f"SV_Insertion\t{sv_ins_count}\t{rate(sv_ins_count, total_block_length):.6e}\t{rate(sv_ins_count, effective_block_length):.6e}\n")
    out_f.write(f"SV_Deletion\t{sv_del_count}\t{rate(sv_del_count, total_block_length):.6e}\t{rate(sv_del_count, effective_block_length):.6e}\n")
    out_f.write(f"Tandem_Repeat\t{trf_count}\t{rate(trf_count, total_block_length):.6e}\t{rate(trf_count, effective_block_length):.6e}\n")
    out_f.write(f"Total\t{total_variants}\t{rate(total_variants, total_block_length):.6e}\t{rate(total_variants, effective_block_length):.6e}\n")
    out_f.write(f"Block_length_bp\t{total_block_length}\t.\t.\n")
    out_f.write(f"Block_length_effective_bp\t{effective_block_length}\t.\t.\n")

# --- WRITE INDEL LENGTHS FOR HISTOGRAM ---
with open(output_lengths_tsv, 'w') as len_f:
    len_f.write("Indel_Length\n")
    for length in indel_lengths:
        len_f.write(f"{length}\n")

print(f"Done. Output written to:")
print(f"  {output_stats_tsv}")
print(f"  {output_lengths_tsv}")
print(f"  {output_bed}")
