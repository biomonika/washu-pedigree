#!/bin/bash

# ============================================================
# USER CONFIGURATION — edit these paths before running
# (chrX is replaced by the actual chromosome at runtime)
# ============================================================
BLOCKS_DIR="/path/to/paternal_transmitted_blocks_per_chr"   # folder with per-chr block TSVs
CACTUS_RESULTS_DIR="/path/to/cactus/results"                # folder containing PAN027.chrN.project_PAN027_pat/ subfolders
OUTPUT_DIR="/path/to/output"                                 # where filtered VCFs will be written
ANNOTATIONS_DIR="/path/to/annotations"                       # folder with BED annotation files

PROBLEMATIC_PAN027="${ANNOTATIONS_DIR}/problematic.PAN027.bed"
PROBLEMATIC_PAN011="${ANNOTATIONS_DIR}/problematic.PAN011.bed"
SW_PRONE="${ANNOTATIONS_DIR}/sw_prone_regions.bed"
# ============================================================

BLOCKS="${BLOCKS_DIR}/chrX.tsv"
INPUT_VCF="${CACTUS_RESULTS_DIR}/PAN027.chrX.project_PAN027_pat/PAN027.chrX.project_PAN027_pat.wave.vcf.gz"
GBZ="${CACTUS_RESULTS_DIR}/PAN027.chrX.project_PAN027_pat/PAN027.chrX.project_PAN027_pat.gbz"
OUTPUT_VCF="${OUTPUT_DIR}/PAN027.chrX.project_PAN027_pat.wave.filtered.vcf"
OUTPUT_VCF_UNFILTERED="${OUTPUT_DIR}/PAN027.chrX.project_PAN027_pat.wave.unfiltered.vcf"
# ============================================================
# To reduce filtering stringency, comment out lines in STEP 6:
#   - Comment out EXCLUDE_B block: skip PAN011 problematic regions
#   - Comment out EXCLUDE_A block: skip PAN027 problematic regions
#   (SW-prone / EXCLUDE_C is always recommended to keep)
# ============================================================

# Automatically detect REF path from GBZ
REF_PATH=$(vg paths -L -x "$GBZ" 2>/dev/null | grep "PAN027_pat" | grep -v "#[0-9]*$" | head -1)
if [[ -z "$REF_PATH" ]]; then
    echo "ERROR: Could not find PAN027_pat reference path in GBZ."
    exit 1
fi
echo "Using reference path: $REF_PATH"

mkdir -p "$(dirname "$OUTPUT_VCF")"

# Path to Python helper (must be in the same directory as this script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_LOOKUP="${SCRIPT_DIR}/vg_pos_lookup.py"

# Temporary files - auto-cleaned on exit
TEMP_VCF=$(mktemp --suffix=.vcf)
TEMP_VCF2=$(mktemp --suffix=.vcf)
TEMP_VCF3=$(mktemp --suffix=.vcf)
POS_INPUT=$(mktemp)
POS_RESULTS=$(mktemp)
BED_FILE=$(mktemp --suffix=.bed)
VG_JSON=$(mktemp --suffix=.json)
BED_REF_HITS=$(mktemp --suffix=.bed)
BED_PAN011_HITS=$(mktemp --suffix=.bed)

trap 'rm -f "$TEMP_VCF" "$TEMP_VCF2" "$TEMP_VCF3" "$POS_INPUT" "$POS_RESULTS" "$BED_FILE" "$VG_JSON" "$BED_REF_HITS" "$BED_PAN011_HITS"' EXIT

# Check required tools
if ! command -v bcftools &> /dev/null; then
    echo "ERROR: bcftools not found in PATH."
    exit 1
fi
if ! command -v vg &> /dev/null; then
    echo "ERROR: vg not found in PATH."
    exit 1
fi
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found in PATH."
    exit 1
fi
if [[ ! -f "$PY_LOOKUP" ]]; then
    echo "ERROR: Python helper not found: $PY_LOOKUP"
    exit 1
fi

# ============================================================
# STEP 1: PREPARE VCF HEADER FOR 2 SAMPLES
# Paternal inheritance line: grandfather (PAN011) -> father (PAN027) -> granddaughter (PAN028)
# ============================================================
bcftools view -h "$INPUT_VCF" | grep -v "#CHROM" > "$TEMP_VCF"

cat >> "$TEMP_VCF" << 'EOF'
##INFO=<ID=POS_PAN028,Number=1,Type=Integer,Description="Position of variant in PAN028 haplotype contig">
##INFO=<ID=POS_PAN011,Number=1,Type=Integer,Description="Position of variant in PAN011 haplotype contig">
##INFO=<ID=CONTIG_PAN028,Number=1,Type=String,Description="PAN028 haplotype contig name">
##INFO=<ID=CONTIG_PAN011,Number=1,Type=String,Description="PAN011 haplotype contig name">
EOF

echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPAN028\tPAN011" >> "$TEMP_VCF"

echo "------------------------------------------------------------"
echo "STEP 1: Filtering variants and extracting 2 sample columns..."
echo "------------------------------------------------------------"

while read -r line; do
    p11_col=$(echo "$line" | awk '{print $1}')
    p27_col=$(echo "$line" | awk '{print $2}')
    p28_col=$(echo "$line" | awk '{print $3}')

    if [[ "$p11_col" == *"haplotype1"* ]]; then p11_idx=3; p11_hap="PAN011_hap1"; else p11_idx=4; p11_hap="PAN011_hap2"; fi
    if [[ "$p28_col" == *"haplotype1"* ]]; then p28_idx=6; p28_hap="PAN028_hap1"; else p28_idx=7; p28_hap="PAN028_hap2"; fi

    region=$(echo "$p27_col" | cut -d':' -f2)
    chrom=$(echo "$p27_col" | cut -d':' -f1)

    # Retain variants where granddaughter (PAN028) is REF and grandfather (PAN011) carries the allele
    bcftools view -H -r "${chrom}:${region}" "$INPUT_VCF" 2>/dev/null | \
    awk -v p28=$((p28_idx + 10)) -v p11=$((p11_idx + 10)) \
        -v p28h="$p28_hap" -v p11h="$p11_hap" \
        -v pos_input="$POS_INPUT" \
        -v vcf_out="$TEMP_VCF" \
    '
    {
        split($p28, a, ":"); gt28 = a[1];
        split($p11, b, ":"); gt11 = b[1];

        if (gt28 == "0" && gt11 ~ /^[1-9]/) {
            print $2 "\t" p28h "\t" p11h >> pos_input;
            for (i=1; i<=9; i++) printf "%s\t", $i >> vcf_out;
            printf "%s\t%s\n", $p28, $p11 >> vcf_out;
        }
    }'

done < "$BLOCKS"

echo "------------------------------------------------------------"
echo "STEP 2: Building BED and running vg find (single call)..."
echo "------------------------------------------------------------"

awk -v ref="$REF_PATH" '{print ref "\t" ($1-1) "\t" $1}' "$POS_INPUT" > "$BED_FILE"

vg find -x "$GBZ" -R "$BED_FILE" 2>/dev/null | vg view -j - 2>/dev/null > "$VG_JSON"

echo "------------------------------------------------------------"
echo "STEP 3: Computing ALT haplotype positions from vg JSON..."
echo "------------------------------------------------------------"

python3 "$PY_LOOKUP" "$VG_JSON" "$POS_INPUT" "$REF_PATH" > "$POS_RESULTS"

echo "------------------------------------------------------------"
echo "STEP 4: Inserting position tags into VCF INFO column..."
echo "------------------------------------------------------------"

grep "^#" "$TEMP_VCF" > "$TEMP_VCF2"

paste <(grep -v "^#" "$TEMP_VCF") "$POS_RESULTS" | \
awk 'BEGIN{OFS="\t"} {
    pos28 = $(NF-3);
    con28 = $(NF-2);
    pos11 = $(NF-1);
    con11 = $NF;

    extra = "";
    if (pos28 != ".") extra = extra ";POS_PAN028=" pos28 ";CONTIG_PAN028=" con28;
    if (pos11 != ".") extra = extra ";POS_PAN011=" pos11 ";CONTIG_PAN011=" con11;
    if (extra != "") $8 = $8 extra;

    for (i=1; i<=NF-4; i++) printf "%s%s", $i, (i==NF-4 ? "\n" : "\t");
}' >> "$TEMP_VCF2"

mv "$TEMP_VCF2" "$TEMP_VCF"

echo "------------------------------------------------------------"
echo "STEP 5: ALT allele trimming and deduplication..."
echo "------------------------------------------------------------"

bcftools view -a "$TEMP_VCF" | bcftools norm -d exact -o "$TEMP_VCF3"

cp "$TEMP_VCF3" "$OUTPUT_VCF_UNFILTERED"
echo "Unfiltered output saved: $OUTPUT_VCF_UNFILTERED"

# ============================================================
# STEP 6: EXCLUDE VARIANTS IN PROBLEMATIC / SW-PRONE REGIONS
# Logic: remove a variant if ANY of the following is true (OR):
#   (A) REF position overlaps problematic.PAN027.bed        → EXCLUDE_A
#   (B) PAN011 alt position overlaps problematic.PAN011.bed → EXCLUDE_B
#   (C) REF position overlaps sw_prone_regions.bed           → EXCLUDE_C
#
# To reduce filtering, comment out EXCLUDE_A and/or EXCLUDE_B blocks below.
# ============================================================
echo "------------------------------------------------------------"
echo "STEP 6: Removing variants in problematic / switch error-prone regions..."
echo "------------------------------------------------------------"

grep "^#" "$TEMP_VCF3" > "$OUTPUT_VCF"

grep -v "^#" "$TEMP_VCF3" | \
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, NR}' > "${BED_REF_HITS}"

grep -v "^#" "$TEMP_VCF3" | \
awk 'BEGIN{OFS="\t"} {
    contig = "."; pos = ".";
    n = split($8, info, ";");
    for (i=1; i<=n; i++) {
        if (info[i] ~ /^CONTIG_PAN011=/) contig = substr(info[i], index(info[i],"=")+1);
        if (info[i] ~ /^POS_PAN011=/)   pos    = substr(info[i], index(info[i],"=")+1);
    }
    if (contig != "." && pos != ".") print contig, pos-1, pos, NR;
}' > "${BED_PAN011_HITS}"

# (A) REF pos overlaps problematic.PAN027.bed
EXCLUDE_A=$(bedtools intersect -a "$BED_REF_HITS" -b "$PROBLEMATIC_PAN027" -wa \
            | awk '{print $4}' | sort -u)

# (B) PAN011 alt pos overlaps problematic.PAN011.bed
EXCLUDE_B=""
if [[ -s "$BED_PAN011_HITS" ]]; then
    EXCLUDE_B=$(bedtools intersect -a "$BED_PAN011_HITS" -b "$PROBLEMATIC_PAN011" -wa \
                | awk '{print $4}' | sort -u)
fi

# (C) REF pos overlaps sw_prone_regions.bed
EXCLUDE_C=$(bedtools intersect -a "$BED_REF_HITS" -b "$SW_PRONE" -wa \
            | awk '{print $4}' | sort -u)

ALL_EXCLUDE=$(printf "%s\n%s\n%s\n" "$EXCLUDE_A" "$EXCLUDE_B" "$EXCLUDE_C" \
              | grep -v '^$' | sort -nu)

TOTAL_BEFORE=$(grep -v "^#" "$TEMP_VCF3" | wc -l)

if [[ -z "$ALL_EXCLUDE" ]]; then
    echo "  No variants overlap any problematic/sw-prone region — nothing removed."
    grep -v "^#" "$TEMP_VCF3" >> "$OUTPUT_VCF"
else
    grep -v "^#" "$TEMP_VCF3" | \
    awk -v excl="$ALL_EXCLUDE" '
    BEGIN {
        n = split(excl, arr, "\n");
        for (i=1; i<=n; i++) bad[arr[i]] = 1;
    }
    { if (!bad[NR]) print }
    ' >> "$OUTPUT_VCF"

    TOTAL_AFTER=$(grep -v "^#" "$OUTPUT_VCF" | wc -l)
    REMOVED=$(( TOTAL_BEFORE - TOTAL_AFTER ))
    echo "  Removed ${REMOVED} of ${TOTAL_BEFORE} variants (PAN027_problematic: $(echo "$EXCLUDE_A" | grep -c .) | PAN011_problematic: $(echo "$EXCLUDE_B" | grep -c .) | SW-prone: $(echo "$EXCLUDE_C" | grep -c .) — some may overlap)."
    echo "  Retained ${TOTAL_AFTER} variants."
fi

echo "------------------------------------------------------------"
echo "DONE!"
echo "  Unfiltered: $OUTPUT_VCF_UNFILTERED"
echo "  Filtered:   $OUTPUT_VCF"
echo "------------------------------------------------------------"
