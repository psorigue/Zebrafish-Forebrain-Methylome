
file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/5mC_all.txt"

ds <- read.csv(file, header = F, sep = "\t")

colnames(ds) <- c("chr", "start", "end", "name", "meth", "cov", "nCpG_cov", "replicate", "type")

library(ggplot2)

p <- ggplot(ds, aes(x = type, y = meth)) +
  geom_violin(trim = T, fill = "steelblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  labs(
    x = "Genomic feature type",
    y = "Methylation level",
    title = "Methylation distribution by genomic type"
  ) +
  theme_classic()


out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/5mC_all.pdf"
ggsave(out_file, p, device = "pdf")
