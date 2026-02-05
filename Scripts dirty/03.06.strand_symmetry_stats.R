#====================================================
# Replicate-aware strand comparison for CH methylation
#====================================================

library(dplyr)
library(tidyr)
library(ggplot2)

#=============================
# 1. Genomic bins
#=============================

name_analysis <- "genome_50kb_bins"
file <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/output/", 
               name_analysis, "/CH_", name_analysis, "_stranded.txt")
ds <- read.csv(file, header = FALSE, sep = "\t")

# Set column names
colnames(ds) <- c("chr", "start", "end", "name", "meth", "cov", "nCpG_cov", "replicate", "strand", "motif")

#------------------------------------
# Step 1: compute per-replicate mean methylation
#------------------------------------
ds_rep <- ds %>%
  group_by(chr, motif, strand, replicate) %>%
  summarise(
    mean_meth = mean(meth, na.rm = TRUE),
    .groups = "drop"
  )

#------------------------------------
# Step 2: t-test per chr × motif across replicates
#------------------------------------
chr_motif_results <- ds_rep %>%
  group_by(chr, motif) %>%
  summarise(
    n_pos = sum(strand == "pos"),
    n_neg = sum(strand == "neg"),
    mean_pos = mean(mean_meth[strand == "pos"], na.rm = TRUE),
    mean_neg = mean(mean_meth[strand == "neg"], na.rm = TRUE),
    mean_diff = mean_pos - mean_neg,
    p_value = ifelse(
      n_pos > 1 & n_neg > 1,
      t.test(mean_meth[strand == "pos"], mean_meth[strand == "neg"], var.equal = FALSE)$p.value,
      NA
    ),
    .groups = "drop"
  ) %>%
  group_by(motif) %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Write output
out_file <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/stats/", 
                   name_analysis, "_strand_stats.txt")
write.table(chr_motif_results, out_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


#=============================
# 2. Promoters and genes
#=============================

name_analysis <- "promoters"
file <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/output/", 
               name_analysis, "/CH_", name_analysis, "_stranded.txt")
ds <- read.csv(file, header = FALSE, sep = "\t")

# Set column names
colnames(ds) <- c("chr", "start", "end", "name", "dum", "st_region", "meth", "cov", 
                  "nCpG_cov", "replicate", "strand", "motif")

#------------------------------------
# Step 1: compute per-replicate mean methylation per motif × strand
#------------------------------------
ds_rep <- ds %>%
  group_by(motif, strand, replicate) %>%
  summarise(
    mean_meth = mean(meth, na.rm = TRUE),
    .groups = "drop"
  )

#------------------------------------
# Step 2: t-test per motif across replicates
#------------------------------------
motif_results_all <- ds_rep %>%
  group_by(motif) %>%
  summarise(
    n_pos = sum(strand == "pos"),
    n_neg = sum(strand == "neg"),
    mean_pos = mean(mean_meth[strand == "pos"], na.rm = TRUE),
    mean_neg = mean(mean_meth[strand == "neg"], na.rm = TRUE),
    mean_diff = mean_pos - mean_neg,
    p_value = ifelse(
      n_pos > 1 & n_neg > 1,
      t.test(mean_meth[strand == "pos"], mean_meth[strand == "neg"], var.equal = FALSE)$p.value,
      NA
    ),
    .groups = "drop"
  ) %>%
  mutate(padj = p.adjust(p_value, method = "BH"))

# Write output
out_file <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/stats/", 
                   name_analysis, "_strand_stats.txt")
write.table(motif_results_all, out_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
