# For XY sample, set chrX genotypes (default half-missing) to homozygous
# so that flare sees chrX as two identical haplotypes

#================#
# Get input data #
#================#

library(optparse)

# Get command line arguments
option_list <- list(
  make_option(c("--input"), type = "character", default = NULL, 
              help = "Input VCF file path", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL, 
              help = "Output VCF file path", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("Please provide input VCF path and output VCF path as command line arguments.")
}

vcfpath <- opt$input
output_vcfpath <- opt$output


#========================#
# Read in and modify VCF #
#========================#

library(tidyverse)
library(data.table)
library(pbmcapply)

# Read in VCF file
vcf <- fread(cmd = paste0("zcat ", vcfpath,
                          " | grep -v '##'"))

homozygose_gts <- function(index, in_table) {
  subset <- in_table[index, ]
  old_col <- subset$syndip
  
  # Get original genotype field
  old_gt <- strsplit(old_col, ":")[[1]][1]
  # If just one of the GTs is missing, convert it to homozygous
  if ((old_gt != "./.") & (old_gt != ".|.") & (old_gt %like% "\\.")) {
    # Get the GT value
    haploid_gt <- gsub("[^0-9]", "", old_gt)
    # Create false homozygous genotype
    new_gt <- paste0(haploid_gt, "|", haploid_gt)
  } else {
    # Otherwise, keep the GT the same
    new_gt <- old_gt
  }
  
  # Add on everything after genotype field if it exists
  other_gt_fields <- strsplit(old_col, ":")[[1]][-1]
  if (length(other_gt_fields > 0)) {
    other_gt_fields <- other_gt_fields %>%
      paste0(., collapse = ":")
    new_col <- paste0(new_gt, ":", other_gt_fields)
  } else {
    new_col <- new_gt
  }
  
  subset[, syndip := new_col]
  return(subset)
}
outvcf <- pbmclapply(1:nrow(vcf),
                     function(x) homozygose_gts(x, vcf),
                     mc.cores = 16) %>%
  rbindlist()


#===============#
# Write new VCF #
#===============#

# Output old VCF header
system(command = paste0("zcat ", vcfpath, " | grep '##'",
                        " > ", output_vcfpath))
# Write VCF
fwrite(outvcf, file = output_vcfpath,
       sep = "\t", col.names = TRUE, append = TRUE)
