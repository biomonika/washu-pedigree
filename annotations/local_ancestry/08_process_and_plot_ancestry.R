# Code to merge flare local ancestry calls into regions
# and plot local ancestry along each chromosome

library(tidyverse)
library(data.table)
library(khroma)
library(wesanderson)
library(GenomicRanges)
library(pbmcapply)

setwd("/scratch16/rmccoy22/syan11/washu_local_ancestry/08_ancestry_plots")

sample <- "PAN028"
# Haplotype names for this sample
hap_names <- c("haplotype1", "haplotype2")

# Path to flare output vcfs
flare_path <- "/scratch16/rmccoy22/syan11/washu_local_ancestry/08_ancestry_plots/flare_lifted/"
# Directory with censat annotations and assembly .fai files (to get chr lengths)
assembly_dir <- "pedigree_assemblies_v1.0_censats/"

# Path to censat annotations
censat_path <- paste0(assembly_dir, "assembly.v1.0.", sample, ".diploid.cenSat.bed")
# Handle only one chrX haplotype for PAN011
hap2_chrs <- c(1:22, "X_par1", "X_nonpar", "X_par2")
if (sample == "PAN011") {
  hap2_chrs <- c(1:22)
}


#====================#
# Read in flare VCFs #
#====================#

read_data <- function(sample, i, haplotype, hap_to_keep) {
  path <- paste0(flare_path, sample, ".chr", i, ".anc.", haplotype, "_lifted.vcf.gz")
  
  hap <- fread(path) %>%
    # split FORMAT column
    .[, c("GT", "AN1", "AN2", "ANP1", "ANP2") := tstrsplit(syndip, split = ":")] %>%
    # combine info for haplotype
    .[, hap1 := paste0(AN1, ",", ANP1)] %>%
    .[, hap2 := paste0(AN2, ",", ANP2)] %>%
    .[, -c("ID", "QUAL", "FILTER", "INFO", "FORMAT", "syndip", "GT", "AN1", "AN2", "ANP1", "ANP2")] %>%
    # pivot so that haplotypes are separate rows
    pivot_longer(cols = c(hap1, hap2),
                 names_to = "haplotype", values_to = "tmp") %>%
    setDT() %>%
    # make separate columns for ancestry and probabilities of AFR/EUR
    .[, c("ancestry", "p_EUR", "p_AFR") := tstrsplit(tmp, ",")] %>%
    .[haplotype == hap_to_keep] %>%
    .[, -c("tmp", "haplotype")] %>%
    # make END column that goes to the position of the next SNP
    group_by(`#CHROM`) %>%
    mutate(END = lead(POS)) %>%
    ungroup() %>%
    setDT() %>%
    # set END position for last SNP in each chromosome
    .[is.na(END), END := POS]
}
hap1 <- pbmclapply(as.list(c(1:22, "X_par1", "X_nonpar", "X_par2")),
                   function(x) read_data(sample, x, hap_names[1], "hap1"),
                   mc.cores = 4) %>%
  rbindlist()
hap2 <- pbmclapply(as.list(hap2_chrs),
                   function(x) read_data(sample, x, hap_names[2], "hap2"),
                   mc.cores = 4) %>%
  rbindlist()

# merge intervals for combinations of haplotypes and ancestry probabilities
group_merge <- function(anc, input, hap_mod) {
  # subset to haplotype and ancestry probability of interest
  dt <- input[p_AFR == anc]
  
  # convert to granges object
  gr <- makeGRangesFromDataFrame(dt,
                                 seqnames.field = c("#CHROM"),
                                 start.field = "POS",
                                 end.field = "END") %>%
    reduce()
  # merge intervals
  reduced <- as.data.table(gr)[, -c("width", "strand")] %>%
    .[, p_AFR := anc] %>%
    # set chromosome numbers for plotting
    .[, chr_num := as.numeric(gsub("chr", "", seqnames))] %>%
    .[seqnames == "chrX", chr_num := 23] %>%
    # Complicated calculation to set y-axis position in final plot
    .[, chr_num := chr_num*2 + hap_mod+chr_num-1]
  return(reduced)
}
out1 <- pbmclapply(as.list(unique(hap1$p_AFR)),
                   function(x) group_merge(x, hap1, 1),
                   mc.cores = 4) %>%
  rbindlist()
out2 <- pbmclapply(as.list(unique(hap2$p_AFR)),
                   function(x) group_merge(x, hap2, 0),
                   mc.cores = 4) %>%
  rbindlist()

# Output ancestry bed files
setorder(out1, seqnames, start)
setorder(out2, seqnames, start)
fwrite(out1,
       paste0(sample, ".ancestry.", hap_names[1], ".assemblyCoords.bed"),
       sep = "\t")
fwrite(out2,
       paste0(sample, ".ancestry.", hap_names[2], ".assemblyCoords.bed"),
       sep = "\t")

# Calculate ancestry proportions for each haplotype
lapply(list(out1, out2),
       function(dt) {
         dt %>%
           .[, length := end - start] %>%
           .[, bases_AFR := as.numeric(p_AFR) * length] %>%
           .[, bases_EUR := (1-as.numeric(p_AFR)) * length] %>%
           summarize(prop_AFR = sum(bases_AFR) / sum(length),
                     prop_EUR = sum(bases_EUR) / sum(length))
       })


#=============================#
# Read in and format metadata #
#=============================#

# chromosome lengths (from assembly faidx)
chrom_lengths <- lapply(as.list(hap_names),
                        function(x) fread(paste0(assembly_dir,
                                                 "assembly.v1.0.", sample, ".", x, ".fa.fai")) %>%
                          .[, hap := x]) %>%
  rbindlist() %>%
  .[, c("V1", "V2", "hap")] %>%
  # Get chromosome number from `chr` column
  .[, y_pos := tstrsplit(V1, "r", keep = 2)] %>%
  .[, y_pos := as.numeric(y_pos)] %>%
  .[V1 == "chrX", y_pos := 23] %>%
  # Complicated calculation to set y-axis position in final plot
  .[hap == hap_names[1], y_pos := (y_pos*2 + y_pos+1-1)] %>%
  .[hap == hap_names[2], y_pos := (y_pos*2 + y_pos+0-1)]

# Read in centromere annotations
centromeres <- fread(censat_path,
                     skip = 1) %>%
  # Get chromosome number from `chr` column
  .[, chr := tstrsplit(V1, split = "\\.", keep = 2)] %>%
  .[, y_pos := as.numeric(gsub("chr", "", chr))] %>%
  .[chr == "chrX", y_pos := 23] %>%
  # Complicated calculation to set y-axis position in final plot
  .[V1 %like% hap_names[1], y_pos := (y_pos*2 + y_pos+1-1)] %>%
  .[V1 %like% hap_names[2], y_pos := (y_pos*2 + y_pos+0-1)]


#================#
# Plot ideograms #
#================#

# width of one chromosome bar
width <- 0.75
# chromosome ticks
ticks <- chrom_lengths %>%
  group_by(V1) %>%
  summarize(tick = mean(y_pos)) %>%
  dplyr::pull(., tick) %>%
  sort()
# data in order of appearance on "y" axis of plot
order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
           "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
           "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
           "chr22", "chrX")

ggplot() +
  # ancestry
  geom_rect(data = out1,
            aes(xmin = as.numeric(chr_num) - width/2,
                xmax = as.numeric(chr_num) + width/2,
                ymin = start, ymax = end,
                fill = as.numeric(p_AFR))) +
  geom_rect(data = out2,
            aes(xmin = as.numeric(chr_num) - width/2,
                xmax = as.numeric(chr_num) + width/2,
                ymin = start, ymax = end,
                fill = as.numeric(p_AFR))) +
  scale_fill_gradient2(low = wes_palette("Royal2")[3],
                       mid = wes_palette("Royal2")[4],
                       high = wes_palette("Royal2")[5],
                       midpoint = 0.5) +
  labs(y = "Position (bp)", x = "Chromosome", fill = "P(AFR ancestry)") +
  
  # centromeres
  geom_rect(data = centromeres,
            aes(xmin = as.numeric(y_pos) - width/2,
                xmax = as.numeric(y_pos) + width/2,
                ymin = V2, ymax = V3),
            fill = "gray85") +
  # outlined chromosome bar
  geom_rect(data = chrom_lengths,
            aes(xmin = as.numeric(y_pos) - width/2,
                xmax = as.numeric(y_pos) + width/2,
                ymin = 0, ymax = V2),
            fill = NA, color = "black", linewidth = 0.4) +

  ### plot reformatting
  # flip x and y axes
  coord_flip() +
  theme_classic() +
  # change "y" axis ticks to discrete scale so each chromosome gets its own tick
  scale_x_continuous(breaks = ticks,
                     labels = order,
                     # spacing between "y" axis and start of chromosome bars
                     expand = c(0, width*1.5)) +
  # spacing between "x" axis and bottom of chr1 bar
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Position (bp)", x = "Chromosome", fill = "P(AFR ancestry)") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))
ggsave(paste0("ancestry_", sample, "_assemblyCoords.pdf"),
       width = 11, height = 7.5)