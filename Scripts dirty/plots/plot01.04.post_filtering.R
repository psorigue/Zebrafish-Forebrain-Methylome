library(ggplot2)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1D.post_filter/dataset.txt"
out_pdf <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1D.post_filter/coverage_chr.pdf"
out_tiff <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1D.post_filter/coverage_chr.tiff"

ds <- read.csv(file, header = T, sep = "\t")


# Exclude chromosome NC_002333.2
X_filtered <- ds %>%
  filter(X.rname != "NC_002333.2")

# Replace chromosome names
mapping <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Epi/Ref_genome/GRCz12tu/chr_mapping.tsv",
                    header = F, col.names = c("acc", "name"), sep = "\t")
X_filtered_map <- X_filtered %>%
  left_join(mapping, by = c("X.rname" = "acc"))

# Boxplot of mean depth per chromosome
p <- ggplot(
  X_filtered_map,
  aes(x = reorder(name, meandepth, FUN = median, decreasing = TRUE),
      y = meandepth)
) +
  geom_boxplot(fill = "grey80", width = 0.5) +
  theme_minimal() +
  labs(
    x = "Chromosome",
    y = "Mean Depth",
    title = "Mean Sequencing Depth per Chromosome"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p
ggsave(out_pdf, p)


# Make tiff
{
  tiff(
    out_tiff,
    width = 6,
    height = 4,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p)
  dev.off()
  
}



