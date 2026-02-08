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

ds_region$type <- gsub(ds_region$type, pattern = "50kbp_bin", replacement = "50kbp bin")
ds_region$type <- gsub(ds_region$type, pattern = "CGisland", replacement = "CG island")
ds_region$type <- gsub(ds_region$type, pattern = "gene_body", replacement = "gene body")

mod_colors <- c(
  "5mC"  = "#1f4e99",   # dark blue
  "5hmC" = "#d4a017"    # golden
)

# Plot
# Ensure 'modification' is a factor with the correct order
ds_region$modification <- factor(ds_region$modification, levels = c("5mC", "5hmC"))

dodge_width <- 0.6

p <- ggplot(ds_region, aes(x = type, y = mean_meth, fill = modification)) +
  geom_violin(
    trim = TRUE,
    scale = "area",
    alpha = 0.7,
    position = position_dodge(width = dodge_width)
  ) +
  scale_fill_manual(values = mod_colors) +
  labs(
    x = "Genomic feature",
    y = "Region mean methylation",
    title = "CpG region-level methylation",
    fill = "Modification"
  ) +
  theme_classic()

p

# Make tiff
out_tiff <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/cpg.tiff"

{
  tiff(
    out_tiff,
    width = 10,
    height = 4,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p)
  dev.off()
  
}
