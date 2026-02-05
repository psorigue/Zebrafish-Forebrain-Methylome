library(ggplot2)
library(dplyr)


file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/output/cpg_all.txt"

ds <- read.csv(file, header = T, sep = "\t")


# Summarise data
ds_region <- ds %>%
  group_by(type, modification, name) %>%   # group by region
  summarise(
    mean_meth = mean(meth, na.rm = TRUE), # average across replicates
    .groups = "drop"
  )


# Plot
# Ensure 'modification' is a factor with the correct order
ds$modification <- factor(ds$modification, levels = c("5mC", "5hmC"))

dodge_width <- 0.6

ggplot(ds_region, aes(x = type, y = mean_meth, fill = modification)) +
  geom_violin(
    trim = TRUE,
    scale = "area",
    alpha = 0.7,
    position = position_dodge(width = dodge_width)
  ) +
  geom_boxplot(
    width = 0.04,
    outlier.shape = NA,
    alpha = 0.5,
    position = position_dodge(width = dodge_width)
  ) +
  labs(
    x = "Genomic feature type",
    y = "Mean methylation per region",
    title = "Region-level methylation averaged across replicates",
    fill = "Modification"
  ) +
  theme_classic()


# Statistics
region_stats <- ds_region %>%
  group_by(type, modification) %>%
  summarise(
    median_meth = median(mean_meth, na.rm = TRUE),
    IQR_low     = quantile(mean_meth, 0.25, na.rm = TRUE),
    IQR_high    = quantile(mean_meth, 0.75, na.rm = TRUE),
    IQR         = IQR(mean_meth, na.rm = TRUE),
    n_regions   = n(),
    .groups = "drop"
  )
region_stats

out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/stats/cpg_stats.txt"
write.table(region_stats, file = out_file, col.names = T, row.names = F, sep = "\t", quote = F)
