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

ds_region$type <- gsub(ds_region$type, pattern = "50kbp_bin", replacement = "50kbp bin")
ds_region$type <- gsub(ds_region$type, pattern = "gene_body", replacement = "gene body")


# 2. Genomic bins
ds_bin <- ds_region[ds_region$type == "50kbp bin", ]

# Ensure columns are factors in the correct order
ds_bin$type <- factor(ds_bin$type, levels = unique(ds_bin$type))  # adjust if needed
ds_bin$motif <- factor(ds_bin$motif, levels = c("CA", "CT", "CC"))  # change names to your motifs
ds_bin$strand <- factor(ds_bin$strand, levels = c("pos", "neg"))  # adjust if needed

dodge_width <- 0.8
inner_dodge <- 0.85

strand_colors <- c(
  "pos" = "#1f4e99",   # dark blue
  "neg" = "#d4a017"    # gold
)

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
  scale_fill_manual(
    values = c(
      "pos" = "#1f4e99",
      "neg" = "#d4a017"
    )
  ) +
  labs(
    x = "Motif",
    y = "Methylation level",
    title = "Non-CpG cytosine methylation strandedness",
    fill = "Strand"
  ) +
  coord_cartesian(ylim = c(0, 0.04)) +
  theme_classic()
p1

# Make tiff
out_tiff <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/strand.tiff"

{
  tiff(
    out_tiff,
    width = 5,
    height = 5,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p1)
  dev.off()
  
}


# 3. All types, unstranded
ds_unstranded <- ds_region  # ignore strand
ds_unstranded$motif <- factor(ds_unstranded$motif, levels = c("CA", "CT", "CC"))  # change names to your motifs


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
  scale_fill_manual(
    values = c(
      "CA" = "#1f4e99",
      "CT" = "#d4a017",
      "CC" = "grey80"
    )
  ) +
  labs(
    x = "Genomic feature",
    y = "Region mean methylation",
    title = "Non-CpG cytosine region-level methylation",
    fill = "Motif"
  ) +
  theme_classic()
p2

# Make tiff
out_tiff <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2B.region_mods/noncpg.tiff"

{
  tiff(
    out_tiff,
    width = 6,
    height = 5,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p2)
  dev.off()
  
}
