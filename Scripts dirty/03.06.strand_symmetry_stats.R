library(dplyr)
library(tidyr)


# 1. Genomic bins
name_analysis <- "genome_50kb_bins"
file <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/output/", name_analysis, "/CH_", name_analysis, "_stranded.txt")
ds <- read.csv(file, header = F, sep = "\t")

# If genomic bins
colnames(ds) <- c("chr", "start", "end", "name", "meth", "cov", "nCpG_cov", "replicate", "strand", "motif")

# Do the mean methylation level differ for each chromosome?
chr_motif_results <- ds %>%
  group_by(chr, motif) %>%
  summarise(
    n_pos = sum(strand == "pos"),
    n_neg = sum(strand == "neg"),
    mean_pos = mean(meth[strand == "pos"], na.rm = TRUE),
    mean_neg = mean(meth[strand == "neg"], na.rm = TRUE),
    mean_diff = mean_pos - mean_neg,
    p_value = t.test(
      meth[strand == "pos"],
      meth[strand == "neg"],
      var.equal = FALSE
    )$p.value,
    .groups = "drop"
  )

chr_motif_results

# Correct FDR
chr_motif_results <- chr_motif_results %>%
  group_by(motif) %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Write output
out_file <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/", name_analysis, "_strand_stats.txt")
write.table(chr_motif_results, out_file, sep ="\t", col.names = T, row.names = F,quote = F)


# 2. Promoters and genes
name_analysis <- "promoters"
file <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/output/", name_analysis, "/CH_", name_analysis, "_stranded.txt")
ds <- read.csv(file, header = F, sep = "\t")
colnames(ds) <- c("chr", "start", "end", "name", "dum", "st_region", "meth", "cov", "nCpG_cov", "replicate", "strand", "motif")

motif_results_all <- ds %>%
  group_by(motif) %>%
  summarise(
    n_pos = sum(strand == "pos"),
    n_neg = sum(strand == "neg"),
    mean_pos = mean(meth[strand == "pos"], na.rm = TRUE),
    mean_neg = mean(meth[strand == "neg"], na.rm = TRUE),
    mean_diff = mean_pos - mean_neg,
    p_value = t.test(
      meth[strand == "pos"],
      meth[strand == "neg"],
      var.equal = FALSE
    )$p.value,
    .groups = "drop"
  )

motif_results_all

# Correct FDR
motif_results_all <- motif_results_all %>%
  group_by(motif) %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

out_file <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/", name_analysis, "_strand_stats.txt")
write.table(motif_results_all, out_file, sep ="\t", col.names = T, row.names = F,quote = F)
