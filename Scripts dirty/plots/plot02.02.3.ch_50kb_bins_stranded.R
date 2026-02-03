library(ggplot2)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/CH_genome_50kb_bins_stranded.txt"

ds <- read.csv(file, header = F, sep = "\t")

colnames(ds) <- c("chr", "start", "end", "name", "meth", "cov", "nCpG_cov", "replicate", "strand", "motif")


pd <- position_dodge(width = 0.8)
p <- ggplot(ds, aes(x = motif, y = meth, fill = strand)) +
  geom_violin(scale = "width", width = 0.5, position = pd) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.5,
    position = pd
  ) +
  labs(
    x = "nonCpG Motif",
    y = "Methylation level",
    title = "Methylation distribution",
    fill = "Strand"
  ) +
  theme_classic() +
  coord_cartesian(ylim = c(0, 0.05))
p

out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/CH_genome_50kb_bins_stranded.pdf"
ggsave(out_file, p, device = "pdf")
