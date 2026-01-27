#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})

#args <- c("PAN010.chr17.bed","PAN027.chr17.bed","PAN028.chr17.bed")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript compare_gens_overlap_span.R file1.bed file2.bed file3.bed\n", file=stderr())
  quit(status = 1)
}

# --- helper: read BED as GRanges ---
read_bed_as_gr <- function(path) {
  df <- read.table(path, sep="\t", header=FALSE, quote="", comment.char="",
                   fill=TRUE, stringsAsFactors=FALSE)
  if (ncol(df) < 3)
    stop(sprintf("File %s has fewer than 3 columns", path))
  df <- df[, 1:3]
  names(df) <- c("chrom","start","end")
  bad <- is.na(df$start) | is.na(df$end) | (df$end < df$start)
  if (any(bad))
    df <- df[!bad, , drop=FALSE]
  gr <- GRanges(
    seqnames = df$chrom,
    ranges = IRanges(start = as.integer(df$start) + 1L,
                     end   = as.integer(df$end))
  )
  # --- check chromosome count ---
  chroms <- unique(as.character(seqnames(gr)))
  if (length(chroms) != 1L)
    stop(sprintf("Error: %s contains %d chromosomes (%s). This script requires exactly one chromosome per file.",
                 path, length(chroms), paste(chroms, collapse=", ")))
  gr
}


# --- read inputs ---
gr1 <- read_bed_as_gr(args[1])
gr2 <- read_bed_as_gr(args[2])
gr3 <- read_bed_as_gr(args[3])

## ---------- STRICT OVERLAP (actual intervals) ----------
i12  <- intersect(gr1, gr2, ignore.strand=TRUE)
i123 <- intersect(i12, gr3, ignore.strand=TRUE)
strict_overlap_bp <- sum(width(reduce(i123)))

u123 <- reduce(c(gr1, gr2, gr3), ignore.strand=TRUE)
strict_union_bp <- sum(width(u123))

strict_percent <- if (strict_union_bp > 0)
  100 * strict_overlap_bp / strict_union_bp else 0

## ---------- SPAN OVERLAP (fill gaps per file → span per chromosome) ----------

spanify <- function(gr) {
  if (length(gr) == 0L) return(GRanges())
  by_chr <- split(gr, as.character(seqnames(gr)))
  # build one IRanges per chromosome
  seqs <- names(by_chr)
  starts <- vapply(by_chr, function(x) min(start(x)), integer(1))
  ends   <- vapply(by_chr, function(x) max(end(x)),   integer(1))
  GRanges(seqnames = seqs, ranges = IRanges(start = starts, end = ends))
}

sp1 <- spanify(gr1)
sp2 <- spanify(gr2)
sp3 <- spanify(gr3)

# intersection and union of filled spans
si12  <- intersect(sp1, sp2, ignore.strand=TRUE)
si123 <- intersect(si12, sp3, ignore.strand=TRUE)
span_overlap_bp <- sum(width(reduce(si123)))

su123 <- reduce(c(sp1, sp2, sp3), ignore.strand=TRUE)
span_union_bp <- sum(width(su123))

span_percent <- if (span_union_bp > 0)
  100 * span_overlap_bp / span_union_bp else 0

## ---------- ABSOLUTE MIDPOINT DIFFERENCES (spans) ----------
mid <- function(gr) (start(gr) + end(gr)) / 2
mid1 <- mid(sp1)
mid2 <- mid(sp2)
mid3 <- mid(sp3)

midpoint_absolute_distances_abs_delta_gens12_bp <- round(abs(mid2 - mid1))
midpoint_absolute_distances_abs_delta_gens23_bp <- round(abs(mid3 - mid2))

# sub CDRs
total_bp_subCDR_g1 <- sum(width(reduce(gr1, ignore.strand=TRUE)))
total_bp_subCDR_g2 <- sum(width(reduce(gr2, ignore.strand=TRUE)))
total_bp_subCDR_g3 <- sum(width(reduce(gr3, ignore.strand=TRUE)))

total_bp_subCDR_abs_delta_gens12_bp <- abs(total_bp_subCDR_g2 - total_bp_subCDR_g1)
total_bp_subCDR_abs_delta_gens23_bp <- abs(total_bp_subCDR_g3 - total_bp_subCDR_g2)

# CDRs
total_bp_CDR_g1  <- sum(width(sp1))
total_bp_CDR_g2  <- sum(width(sp2))
total_bp_CDR_g3  <- sum(width(sp3))

total_bp_CDR_abs_delta_gens12_bp <- abs(total_bp_CDR_g2 - total_bp_CDR_g1)
total_bp_CDR_abs_delta_gens23_bp <- abs(total_bp_CDR_g3 - total_bp_CDR_g2)

## ---------- OUTPUT ----------
cat(sprintf("strict_subCDR_overlap_percent\t%.2f%%\n", strict_percent))
cat(sprintf("CDR_overlap(span)_percent\t%.2f%%\n",   span_percent))
cat(sprintf("midpoint_absolute_distances_abs_delta_gens12_bp\t%d\n",
            midpoint_absolute_distances_abs_delta_gens12_bp))
cat(sprintf("midpoint_absolute_distances_abs_delta_gens23_bp\t%d\n",
            midpoint_absolute_distances_abs_delta_gens23_bp))
cat(sprintf("total_bp_subCDR_abs_delta_gens12_bp\t%d\n",
            total_bp_subCDR_abs_delta_gens12_bp))
cat(sprintf("total_bp_subCDR_abs_delta_gens23_bp\t%d\n",
            total_bp_subCDR_abs_delta_gens23_bp))
cat(sprintf("total_bp_CDR_abs_delta_gens12_bp\t%d\n",
            total_bp_CDR_abs_delta_gens12_bp))
cat(sprintf("total_bp_CDR_abs_delta_gens23_bp\t%d\n",
            total_bp_CDR_abs_delta_gens23_bp))