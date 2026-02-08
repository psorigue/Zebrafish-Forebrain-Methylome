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
