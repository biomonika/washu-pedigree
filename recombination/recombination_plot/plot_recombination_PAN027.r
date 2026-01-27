suppressMessages({
  library(karyoploteR)
  library(GenomicRanges)
})

base_dir <- "paternal_PAN011/"
outputpdf <- file.path(base_dir, "PAN11_sharedcontigs_PAN027.pdf")

bed_files <- list(
  hap1        = file.path(base_dir, "PAN011_sharedcontigs_to_PAN027-hap1_threegen_merged.bed"),
  hap2        = file.path(base_dir, "PAN011_sharedcontigs_to_PAN027-hap2_threegen_merged.bed"),
  hap1_twogen = file.path(base_dir, "PAN011_sharedcontigs_to_PAN027-hap1_twogen_merged.bed"),
  hap2_twogen = file.path(base_dir, "PAN011_sharedcontigs_to_PAN027-hap2_twogen_merged.bed")
)

markers_file  <- file.path(base_dir, "variantcalling_recombination_paternal.tsv")
phasesets_file <- file.path(base_dir, "recomb_asm-v1.0_paternal_phasesets.tsv")

custom_genome    <- toGRanges("chm13_complete.txt")
custom_cytobands <- toGRanges("chm13v2.0_cytobands_allchrs.txt")

read_bed_as_granges <- function(filepath) {
  if (file.exists(filepath)) {
    df <- read.table(filepath, sep = '\t', stringsAsFactors = FALSE, col.names = c('chr', 'start', 'end', 'contig'), fill = TRUE, fileEncoding = 'latin1')
    if (nrow(df) > 0) {
      return(makeGRangesFromDataFrame(df, seqnames.field = 'chr', start.field = 'start', end.field = 'end'))
    }
  }
  return(NULL)
}

read_marker_table <- function(filepath) {
  read.table(
    filepath,
    sep = "\t",
    stringsAsFactors = FALSE,
    col.names = c("chr", "start", "end", "label"),
    fill = TRUE,
    fileEncoding = "latin1"
  )
}

pdf(outputpdf, width = 10, height = 24)

pp <- getDefaultPlotParams(plot.type = 2)
pp$leftmargin <- 0.1
pp$bottommargin <- 60
pp$data1outmargin <- 28
pp$data1inmargin <- 5
pp$data1height <- 25
pp$data2height <- 25
pp$ideogramheight <- 20
pp$data2inmargin <- 23
pp$data2outmargin <- 1
pp$data3height <- 10
pp$data4height <- 10

kp <- plotKaryotype(
  plot.type = 2,
  genome = custom_genome,
  cytobands = custom_cytobands,
  plot.params = pp
)
kpAddBaseNumbers(kp, tick.dist = 2e7, tick.len = 5, cex = 1.2)
kpDataBackground(kp, data.panel = 1)

#read and plot haplotype data
hap1_ranges        <- read_bed_as_granges(bed_files$hap1)
hap1_twogen_ranges <- read_bed_as_granges(bed_files$hap1_twogen)
hap2_ranges        <- read_bed_as_granges(bed_files$hap2)
hap2_twogen_ranges <- read_bed_as_granges(bed_files$hap2_twogen)

if (!is.null(hap1_ranges)) {
  kpPlotRegions(kp, data = hap1_ranges, col = "#6495ED", r0 = 0, r1 = 1, data.panel = 1)
}
if (!is.null(hap1_twogen_ranges)) {
  kpPlotRegions(kp, data = hap1_twogen_ranges, col = "#6495ED", r0 = 0, r1 = 1, data.panel = 1)
}
if (!is.null(hap2_ranges)) {
  kpPlotRegions(kp, data = hap2_ranges, col = "#9FE2BF", r0 = 0, r1 = 1, data.panel = 1)
}
if (!is.null(hap2_twogen_ranges)) {
  kpPlotRegions(kp, data = hap2_twogen_ranges, col = "#9FE2BF", r0 = 0, r1 = 1, data.panel = 1)
}

# plot markers from variant-based approach
vcpoints <- read_marker_table(markers_file)
vcpoints_ranges <- makeGRangesFromDataFrame(
  vcpoints,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  na.rm = TRUE
)
kpPlotMarkers(
  kp,
  vcpoints_ranges,
  labels = vcpoints$label,
  r0 = 0,
  r1 = 1.3,
  line.color = "#DE3163",
  label.color = "#DE3163",
  text.orientation = "horizontal",
  data.panel = 1
)

# plot phase sets from WhatsHap
ps <- read_marker_table(phasesets_file)
ps_ranges <- makeGRangesFromDataFrame(
  ps,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  na.rm = TRUE
)
kpPlotMarkers(
  kp,
  ps_ranges,
  labels = ps$label,
  r0 = 0,
  r1 = 1.3,
  line.color = "black",
  label.color = "black",
  text.orientation = "horizontal",
  data.panel = 1
)

dev.off()
