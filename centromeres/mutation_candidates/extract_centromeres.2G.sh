#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --job-name="extract_centromeres.2G.20260608"
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --partition=medium
#SBATCH --output=extract_centromeres.20260608.%j.log

set -e
set -x

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/.conda/envs/methylation

TX="transmissions.txt"          # your transmission file (with header chr1 ... chrX)
BED="pedigree.active_hor.bed"  # your BED (col1 == FASTA header)
FA="pedigree.fasta"            # your FASTA

# read chromosome names from header (tab-delimited)
IFS=$'\t' read -r -a CHRS < <(head -n1 "$TX")

# loop each chromosome (i.e., each column)
for idx in "${!CHRS[@]}"; do
  col=$((idx+1))
  chr="${CHRS[$idx]}"

  echo "==> ${chr}"

  # targets for this column (skip header, unique)
  awk -v c="$col" 'NR>1 && $c!=""{print $c}' "$TX" | sort -u > ".targets.${chr}.txt"

  # add tolerance for merging
  # select matching BED rows (col1 in targets) and extract; headers become >TARGET:START-END
  awk -v OFS='\t' 'NR==FNR{t[$0]; next} ($1 in t){print $1,$2,$3,$1}' ".targets.${chr}.txt" "$BED" \
    | bedtools merge -d 500000 | bedtools getfasta -fi "$FA" -bed - -name+ -fo "${chr}.fa"

  rm -f ".targets.${chr}.txt"
done


#split into individual sequences
for f in chr*.fa; do
  [ -e "$f" ] || continue
  base="${f%.fa}"
  echo "==> $f"

  # count headers (= number of sequences) and enforce multiple of 2
  n=$(grep -c '^>' "$f")
  if (( n % 2 )); then
    echo "ERROR: $f has $n sequences (not a multiple of 2)" >&2
  fi

  k=$(( n / 2 ))  # sequences per role

  # Stream through once: first k -> grandparent, next k -> mother
  awk -v base="$base" -v k="$k" '
    # Choose output file based on record index i
    function pick_out(i,   role, idx, suf) {
      if (i <= k)        { role="grandparent";  idx=i }
      else               { role="mother";idx=i-k }
      suf = (k>1) ? "_" idx : ""
      return base "_" suf role ".fa"
    }

    /^>/ {                      # header line -> start a new record
      i++
      out = pick_out(i)
      print $0 > out            # write header to its file
      next
    }
    { print $0 >> out }         # write sequence lines to the same file
  ' "$f"
done
