# ============================================================
# USER CONFIGURATION — edit these paths before running
# ============================================================
CONDA_SH="/opt/miniconda/etc/profile.d/conda.sh"   # path to conda.sh
CONDA_ENV="/private/home/mcechova/conda/ONT"         # conda environment path or name
INPUT_DIR="results_filtered_PAN027_only"             # folder with filtered VCFs (input)
OUTPUT_DIR="results_filtered_PAN027_only_easy_difficult"  # output folder
EASY_BED="PAN027.v1.1.easy.bed"                     # easy regions BED file
DIFFICULT_BED="PAN027.v1.1.difficult.bed"            # difficult regions BED file
# ============================================================

source "$CONDA_SH"
conda activate "$CONDA_ENV"

set -e
#set -x
#take filtered results and split into easy and difficult regions now
#potential overlaps go to difficult regions

for a in "${INPUT_DIR}"/*/*.filtered.vcf; do

    echo "Processing: $a"

    middle_folder=$(basename "$(dirname "$a")")
    b=$(basename "$a")
    vcf="${b%.vcf}"

    outdir="${OUTPUT_DIR}/$middle_folder"
    mkdir -p "$outdir"

    easy_raw="$outdir/${vcf}.easy.raw.vcf.gz"
    diff_vcf="$outdir/${vcf}.difficult.vcf.gz"
    easy_final="$outdir/${vcf}.easy.vcf.gz"

    if [ -f "$easy_raw" ] && [ -f "$diff_vcf" ] && [ -f "$easy_final" ]; then
        echo "Skipping (already exists): $outdir"
        continue
    fi

    # 1. raw splits
    bedtools intersect -header -a "$a" -b "$EASY_BED"      | bcftools sort | bgzip -c > "$easy_raw"
    bedtools intersect -header -a "$a" -b "$DIFFICULT_BED" | bcftools sort | bgzip -c > "$diff_vcf"

    #index the raw bed files
    tabix -p vcf ${easy_raw}
    tabix -p vcf ${diff_vcf}

    # 2. find overlap (variants in both sets)
    easy_raw_count=$(bcftools view -H "$easy_raw" | wc -l)
    diff_count=$(bcftools view -H "$diff_vcf" | wc -l)

    overlap_file=""
    if [ "$easy_raw_count" -gt 0 ] && [ "$diff_count" -gt 0 ]; then
        tmpdir=$(mktemp -d)
        echo "tmpdir: $tmpdir"
        bcftools isec -p "$tmpdir" -Oz \
            "$easy_raw" "$diff_vcf"
        overlap_file="$tmpdir/0002.vcf.gz"
    fi

    # 3. subtract overlap from easy (THIS IS THE KEY STEP)
    overlap_count=0
    if [ -n "$overlap_file" ] && [ -f "$overlap_file" ]; then
        overlap_count=$(bcftools view -H "$overlap_file" | wc -l)
    fi

    if [ "$overlap_count" -gt 0 ]; then
        bcftools isec -C \
            "$easy_raw" "$overlap_file" -w1 \
            -Oz -o "$easy_final"
    else
        cp "$easy_raw" "$easy_final"
    fi

    if [ -n "$tmpdir" ]; then
        rm -rf "$tmpdir"
        tmpdir=""
    fi

    # 4. counts
    orig=$(bcftools view -H "$a" | wc -l)
    easy_before=$(bcftools view -H "$easy_raw" | wc -l)
    easy_after=$(bcftools view -H "$easy_final" | wc -l)
    diff=$(bcftools view -H "$diff_vcf" | wc -l)

    total=$((easy_after + diff))

    echo "  original : $orig"
    echo "  easy     : $easy_before (before removing overlap)"
    echo "  easy     : $easy_after (after removing overlap)"
    echo "  difficult: $diff (includes overlap)"
    echo "  total    : $total"

    if [ "$total" -ne "$orig" ]; then
        echo "WARNING: mismatch in $a"
    fi

    echo "---------------------------"

done