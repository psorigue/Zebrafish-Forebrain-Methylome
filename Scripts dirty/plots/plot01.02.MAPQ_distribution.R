library(tidyverse)
library(scales)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1B.MAPQ_distribution/mapq_distribution.tsv"
out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1B.MAPQ_distribution/mapq_distribution.pdf"

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
  geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = blue_reps) +
  labs(
    x = "MAPQ",
    y = NULL,
    title = "Mapping quality distribution"
  ) +
  theme_classic() +
  theme(legend.position = "none")

pC

ggsave(out_file, pC)