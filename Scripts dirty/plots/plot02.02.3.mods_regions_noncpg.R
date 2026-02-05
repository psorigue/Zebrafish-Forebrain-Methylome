library(ggplot2)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/output/noncpg_all.txt"

ds <- read.csv(file, header = T, sep = "\t")


# Ensure columns are factors in the correct order
ds$type <- factor(ds$type, levels = unique(ds$type))  # adjust if needed
ds$motif <- factor(ds$motif, levels = c("CA", "CC", "CT"))  # change names to your motifs
ds$strand <- factor(ds$strand, levels = c("pos", "neg"))  # adjust if needed

# Create a combined factor for motif + strand
ds$motif_strand <- interaction(ds$motif, ds$strand, sep = "_")

# Plot
dodge_width <- 0.8  # big dodge for type separation
inner_dodge <- 0.7  # inner dodge for motif+strand

ggplot(ds, aes(x = type, y = meth, fill = strand)) +
  geom_violin(
    aes(group = motif_strand), 
    trim = TRUE, 
    alpha = 0.7, 
    position = position_dodge(width = inner_dodge)
  ) +
  geom_boxplot(
    aes(group = motif_strand), 
    width = 0.05, 
    outlier.shape = NA, 
    alpha = 0.5, 
    position = position_dodge(width = inner_dodge)
  ) +
  labs(
    x = "Genomic feature type",
    y = "Methylation level",
    title = "Methylation by type, motif, and strand",
    fill = "Strand"
  ) +
  theme_classic()


out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/CH_genome_50kb_bins_stranded.pdf"
ggsave(out_file, p, device = "pdf")
