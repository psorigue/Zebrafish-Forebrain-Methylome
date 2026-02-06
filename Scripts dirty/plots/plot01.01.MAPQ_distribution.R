library(tidyverse)
library(scales)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/QC/qc_stats/mapq_distribution.tsv"
out_pdf <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1A.MAPQ_distribution/mapq_distribution.pdf"
out_tiff <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1A.MAPQ_distribution/mapq_distribution.tiff"

mapq <- read.csv(file, col.names = c("sample", "mapq"), sep = "\t")

# Ensure sample order is meaningful
mapq$sample <- factor(mapq$sample, levels = sort(unique(mapq$sample)))

# Sober blue palette (light → dark)
blue_reps <- c(
  "#dce9f5",
  "#bcd7ee",
  "#8fbde6",
  "#5fa2da",
  "#2f78c4",
  "#1f4e99"
)

pC <- ggplot(mapq, aes(x = mapq, y = sample, fill = sample)) +
  geom_vline(xintercept = 10, color = "gold", linetype = "dashed", linewidth = 1) +
  geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = blue_reps) +
  labs(
    x = "MAPQ",
    y = NULL,
    title = "Mapping quality distribution"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank()
  )

pC

ggsave(out_pdf, pC)


# Make tiff
{
  tiff(
    out_tiff,
    width = 4,
    height = 4,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  
  print(pC)
  
  dev.off()
  
}
