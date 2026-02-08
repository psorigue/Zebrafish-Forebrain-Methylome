library(ggplot2)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/output/noncpg_all.txt"

ds <- read.csv(file, header = T, sep = "\t")


# 1. Summarise
# Average methylation per region across replicates
ds_region <- ds %>%
  group_by(type, motif, strand, name) %>%  # 'name' = region ID
  summarise(
    mean_meth = mean(meth, na.rm = TRUE),  # average across replicates
    .groups = "drop"
  )

# 2. Genomic bins
ds_bin <- ds_region[ds_region$type == "50kbp_bin", ]

# Stats
violin_stats_bin <- ds_bin %>%
  group_by(motif, strand) %>%
  summarise(
    median_meth = median(mean_meth, na.rm = TRUE),
    IQR_low     = quantile(mean_meth, 0.25, na.rm = TRUE),
    IQR_high    = quantile(mean_meth, 0.75, na.rm = TRUE),
    IQR         = IQR(mean_meth, na.rm = TRUE),
    n_regions   = n(),
    .groups = "drop"
  )
violin_stats_bin
out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/stats/noncpg_genomic_bins_stranded_stats.txt"
write.table(violin_stats_bin, file = out_file, col.names = T, row.names = F, sep = "\t", quote = F)


# 3. All types, unstranded
ds_unstranded <- ds_region  # ignore strand


# Stats
violin_stats_type_motif <- ds_unstranded %>%
  group_by(type, motif) %>%
  summarise(
    median_meth = median(mean_meth, na.rm = TRUE),
    IQR_low     = quantile(mean_meth, 0.25, na.rm = TRUE),
    IQR_high    = quantile(mean_meth, 0.75, na.rm = TRUE),
    IQR         = IQR(mean_meth, na.rm = TRUE),
    n_regions   = n(),
    .groups = "drop"
  )
violin_stats_type_motif
out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/stats/noncpg_unstranded_stats.txt"
write.table(violin_stats_type_motif, file = out_file, col.names = T, row.names = F, sep = "\t", quote = F)



