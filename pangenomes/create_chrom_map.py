import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ============================================================
# USER CONFIGURATION
# ============================================================
BLOCKS_DIR_MAT = "/path/to/maternal_transmitted_blocks_per_chr"
BLOCKS_DIR_PAT = "/path/to/paternal_transmitted_blocks_per_chr"
# ============================================================

# --- ARGUMENTS ---
# Usage: python3 create_chrom_map.py <chr> <hap> [out_base]
if len(sys.argv) >= 3:
    _chr, _hap = sys.argv[1], sys.argv[2]
else:
    _chr, _hap = "1", "pat"
_out_base = sys.argv[3] if len(sys.argv) >= 4 else "results_filtered_all"

_hap_word = "maternal" if _hap == "mat" else "paternal"
_ref      = "PAN027_mat" if _hap == "mat" else "PAN027_pat"

blocks_dir    = BLOCKS_DIR_MAT if _hap == "mat" else BLOCKS_DIR_PAT
blocks_file   = f"{blocks_dir}/chr{_chr}.tsv"
mutations_file = f"{_out_base}/chr{_chr}_{_hap}/heatmap_positions.bed"
output_stacked = f"{_out_base}/chr{_chr}_{_hap}/chromosome_map_stacked.png"
output_alpha   = f"{_out_base}/chr{_chr}_{_hap}/chromosome_map_alpha.png"

# Read chromosome length from VCF header (##contig line)
_vcf_for_len = f"{_out_base}/chr{_chr}_{_hap}/PAN027.chr{_chr}.project_{_ref}.wave.filtered.vcf"
chr_length = None
with open(_vcf_for_len) as _f:
    for _line in _f:
        if not _line.startswith('#'):
            break
        if _line.startswith('##contig') and f'chr{_chr}.{_hap_word}' in _line:
            _m = re.search(r'length=(\d+)', _line)
            if _m:
                chr_length = int(_m.group(1))
                break
if chr_length is None:
    raise ValueError(f"Could not find chromosome length in {_vcf_for_len}")

# --- LOAD BLOCKS ---
blocks = []
with open(blocks_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if len(cols) < 2:
            continue
        match = re.search(r':(\d+)-(\d+)$', cols[1])
        if match:
            start = int(match.group(1))
            end   = int(match.group(2))
            blocks.append((start, end - start))

# --- LOAD MUTATIONS ---
color_map = {
    'SNV':       '#377eb8',  # blue
    'INDEL_1bp': '#984ea3',  # purple
    'INDEL':     '#4daf4a',  # green
    'SV':        '#ff7f00',  # orange
    'TRF':       '#a65628',  # brown
}
track_order = ['SNV', 'INDEL_1bp', 'INDEL', 'SV', 'TRF']

df = None
if os.path.exists(mutations_file):
    df = pd.read_csv(mutations_file, sep='\t', header=None, names=['chr', 'start', 'end', 'type'])
    remap = {
        'INS_1bp': 'INDEL_1bp', 'DEL_1bp': 'INDEL_1bp',
        'INS': 'INDEL', 'DEL': 'INDEL',
        'SV_INS': 'SV', 'SV_DEL': 'SV',
        'MNP': 'SNV',
    }
    df['type'] = df['type'].map(lambda t: remap.get(t, t))


def draw_chromosome_body(ax, y_bottom, track_h):
    ax.broken_barh([(0, chr_length)], (y_bottom, track_h),
                   facecolors='#eaeaea', edgecolors='#bcbcbc')
    if blocks:
        ax.broken_barh(blocks, (y_bottom, track_h),
                       facecolors='#aaaaaa', edgecolors='none', alpha=0.9)


formatter = ticker.FuncFormatter(lambda x, pos: f'{int(x/1e6)}M')

# ============================================================
# VERSION 1: STACKED TRACKS
# ============================================================
n_tracks = len(track_order) + 1  # +1 for chromosome track
track_h  = 0.8
gap      = 0.3
fig_h    = n_tracks * (track_h + gap) + 1.5

fig1, ax1 = plt.subplots(figsize=(14, fig_h))

y_chrom = (n_tracks - 1) * (track_h + gap)
draw_chromosome_body(ax1, y_chrom, track_h)
ax1.text(-chr_length * 0.01, y_chrom + track_h / 2, 'Analyzed blocks',
         va='center', ha='right', fontsize=12, color='#555555')

legend_handles = [plt.Rectangle((0, 0), 1, 1, fc='#aaaaaa', label=f'Analyzed blocks (n={len(blocks)})')]

for i, mut_type in enumerate(track_order):
    y = (n_tracks - 2 - i) * (track_h + gap)
    ax1.broken_barh([(0, chr_length)], (y, track_h),
                    facecolors='#f8f8f8', edgecolors='#dddddd', linewidth=0.5)
    color = color_map.get(mut_type, 'black')
    count = 0
    if df is not None:
        group = df[df['type'] == mut_type]
        count = len(group)
        if count > 0:
            ax1.vlines(group['start'], ymin=y, ymax=y + track_h,
                       colors=color, linewidth=1.0)
    legend_handles.append(plt.Line2D([0], [0], color=color, label=mut_type))
    ax1.text(-chr_length * 0.01, y + track_h / 2, f'{mut_type} (n={count})',
             va='center', ha='right', fontsize=12)

total_h = n_tracks * (track_h + gap)
ax1.set_xlim(0, chr_length)
ax1.set_ylim(-0.5, total_h + 0.5)
ax1.set_yticks([])
ax1.set_xlabel('Chromosome position (bp)', fontsize=17)
ax1.set_title(f'PAN027.chr{_chr}.{_hap_word} – De novo mutations (stacked)', fontsize=17)
ax1.xaxis.set_major_formatter(formatter)
ax1.tick_params(axis='x', labelsize=15)
ax1.legend().set_visible(False)

plt.tight_layout()
plt.savefig(output_stacked, dpi=300, bbox_inches='tight')
plt.close(fig1)
print(f"Stacked: {output_stacked}")

# ============================================================
# VERSION 2: ALPHA OVERLAY (all types on one track)
# ============================================================
fig2, ax2 = plt.subplots(figsize=(14, 2.5))

draw_chromosome_body(ax2, 0, 1)

legend_handles2 = [plt.Rectangle((0, 0), 1, 1, fc='#aaaaaa', label=f'Analyzed blocks (n={len(blocks)})')]

if df is not None:
    for mut_type in track_order:
        color = color_map.get(mut_type, 'black')
        group = df[df['type'] == mut_type]
        if len(group) > 0:
            ax2.vlines(group['start'], ymin=0, ymax=1,
                       colors=color, linewidth=1.2, alpha=0.2)
            legend_handles2.append(plt.Line2D([0], [0], color=color, alpha=0.7,
                                              label=f'{mut_type} (n={len(group)})'))

ax2.set_xlim(0, chr_length)
ax2.set_ylim(-0.5, 1.8)
ax2.set_yticks([])
ax2.set_xlabel('Chromosome position (bp)')
ax2.set_title(f'PAN027.chr{_chr}.{_hap_word} – De novo mutations (alpha overlay)')
ax2.xaxis.set_major_formatter(formatter)
ax2.legend(handles=legend_handles2, loc='upper center',
           bbox_to_anchor=(0.5, -0.35), ncol=len(legend_handles2), frameon=False, fontsize=12)

plt.tight_layout()
plt.savefig(output_alpha, dpi=300, bbox_inches='tight')
plt.close(fig2)
print(f"Alpha: {output_alpha}")
