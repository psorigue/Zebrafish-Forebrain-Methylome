
file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/CH_genome_50kb_bins.txt"

ds <- read.csv(file, header = F, sep = "\t")

colnames(ds) <- c("chr", "start", "end", "name", "meth", "cov", "nCpG_cov", "replicate", "motif")

ggplot(ds, aes(x = motif, y = meth)) +
  geom_violin(trim = T, fill = "steelblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  labs(
    x = "nonCpG Motif",
    y = "Methylation level",
    title = "Methylation distribution"
  ) +
  theme_classic()


p <- ggplot(ds, aes(x = motif, y = meth)) +
  geom_violin(trim = TRUE, fill = "steelblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  scale_y_continuous(limits = c(0, 0.05)) +
  labs(
    x = "nonCpG Motif",
    y = "Methylation level",
    title = "Methylation distribution"
  ) +
  theme_classic()

out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/CH_genome_50kb_bins.pdf"
ggsave(out_file, p, device = "pdf")
