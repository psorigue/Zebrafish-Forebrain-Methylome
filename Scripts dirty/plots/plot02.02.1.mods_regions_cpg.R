library(ggplot2)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/output/cpg_all.txt"

ds <- read.csv(file, header = T, sep = "\t")

# Ensure 'modification' is a factor with the correct order
ds$modification <- factor(ds$modification, levels = c("5mC", "5hmC"))
# Use the same dodge width for violins and boxplots
dodge_width <- 0.6

p <- ggplot(ds, aes(x = type, y = meth, fill = modification)) +
  geom_violin(trim = TRUE, alpha = 0.7, position = position_dodge(width = dodge_width)) +
  geom_boxplot(width = 0.04, outlier.shape = NA, alpha = 0.5, position = position_dodge(width = dodge_width)) +
  labs(
    x = "Genomic feature type",
    y = "Methylation level",
    title = "Methylation distribution by genomic type and modification",
    fill = "Modification"
  ) +
  theme_classic()
p
