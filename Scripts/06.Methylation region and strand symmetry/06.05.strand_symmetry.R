# This script analyzes strand symmetry of methylation in different genomic contexts (genome-wide bins, promoters, genes) by performing t-tests between positive and negative strands across replicates. It also computes mean methylation per strand and motif, and adjusts p-values for multiple testing.

library(dplyr)
library(tidyr)
library(ggplot2)

#=============================
# 1. Genomic bins
#=============================
home <- path.expand("~")
file <- paste0(home, "/Pol/Methylome/methylation_regions/output/genome_50kb_bins/CH_genome_50kb_bins_stranded.txt")
ds <- read.csv(file, header = FALSE, sep = "\t")

# Set column names
colnames(ds) <- c("chr", "start", "end", "name", "meth", "cov", "nCpG_cov", "replicate", "strand", "motif")

# =============================
# 1.1. Chromosome-specific analysis
# =============================

#------------------------------------
# 1.1.1.: Compute per-replicate mean methylation per chr × motif × strand
#------------------------------------
ds_chr_rep <- ds %>%
  group_by(chr, motif, strand, replicate) %>%
  summarise(
    mean_meth = mean(meth, na.rm = TRUE),
    .groups = "drop"
  )

#------------------------------------
# 1.1.2.: Perform t-test per chr × motif 
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
out_file_genomic <- paste0(home, "/Pol/Methylome/methylation_regions/stats/chr_genome_50kb_bins_strand_stats.txt")
write.table(chr_motif_results, out_file_genomic, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# =============================
# 1.2. Genome-wide analysis (no chromosome specificity)
# =============================

#------------------------------------
# 1.2.1.: Compute per-replicate genome-wide mean per motif × strand
#------------------------------------
ds_global_rep <- ds %>%
  group_by(motif, strand, replicate) %>%
  summarise(
    mean_meth = mean(meth, na.rm = TRUE),
    .groups = "drop"
  )

#------------------------------------
# 1.2.2.: Perform t-test per motif
#------------------------------------
global_results <- ds_global_rep %>%
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
  group_by(motif) %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Write output
out_file_genomic <- paste0(home, "/Pol/Methylome/methylation_regions/stats/genome_50kb_bins_strand_stats.txt")
write.table(global_results, out_file_genomic, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


#=============================
# 2. Promoters and genes
#=============================

name_analysis <- "promoters" # promoters, genes 
file <- paste0(home, "/Pol/Methylome/methylation_regions/output/", 
               name_analysis, "/CH_", name_analysis, "_stranded.txt")
ds <- read.csv(file, header = FALSE, sep = "\t")

# Set column names
colnames(ds) <- c("chr", "start", "end", "name", "dum", "st_region", "meth", "cov", 
                  "nCpG_cov", "replicate", "strand", "motif")

#------------------------------------
# 2.1: compute per-replicate mean methylation per motif × strand
#------------------------------------
ds_rep <- ds %>%
  group_by(motif, strand, replicate) %>%
  summarise(
    mean_meth = mean(meth, na.rm = TRUE),
    .groups = "drop"
  )

#------------------------------------
# 2.2: t-test per motif across replicates
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
out_file_genes <- paste0(home, "/Pol/Methylome/methylation_regions/stats/genomic_50kbp_bins_strand_stats.txt")
write.table(motif_results_all, out_file_genes, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
