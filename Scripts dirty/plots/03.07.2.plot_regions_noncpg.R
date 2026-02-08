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

# Ensure columns are factors in the correct order
ds$type <- factor(ds$type, levels = unique(ds$type))  # adjust if needed
ds$motif <- factor(ds$motif, levels = c("CA", "CC", "CT"))  # change names to your motifs
ds$strand <- factor(ds$strand, levels = c("pos", "neg"))  # adjust if needed

dodge_width <- 0.8
inner_dodge <- 0.7

p1 <- ggplot(ds_bin, aes(x = motif, y = mean_meth, fill = strand)) +
  geom_violin(
    trim = TRUE,
    scale = "area",
    alpha = 0.7,
    position = position_dodge(width = inner_dodge)
  ) +
  geom_boxplot(
    width = 0.03,
    outlier.shape = NA,
    alpha = 0.5,
    position = position_dodge(width = inner_dodge)
  ) +
  labs(
    x = "Motif",
    y = "Region-level methylation",
    title = "Methylation by motif and strand (50 kbp bins)",
    fill = "Strand"
  ) +
  coord_cartesian(ylim = c(0, 0.04)) +
  theme_classic()
p1
# Save plot
ggsave("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/CH_genome_50kb_bins_stranded.pdf",
       plot = p1, device = "pdf")


# 3. All types, unstranded
ds_unstranded <- ds_region  # ignore strand

p2 <- ggplot(ds_unstranded, aes(x = type, y = mean_meth, fill = motif)) +
  geom_violin(
    trim = TRUE,
    alpha = 0.7,
    position = position_dodge(width = inner_dodge)
  ) +
  geom_boxplot(
    width = 0.04,
    outlier.shape = NA,
    alpha = 0.5,
    position = position_dodge(width = inner_dodge)
  ) +
  coord_cartesian(ylim = c(0, 0.02)) +
  labs(
    x = "Genomic feature type",
    y = "Region-level methylation",
    title = "Region-level methylation by type and motif",
    fill = "Motif"
  ) +
  theme_classic()
p2
