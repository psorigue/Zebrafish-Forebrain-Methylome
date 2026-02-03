library(tidyverse)
library(scales)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1B.MAPQ_distribution/mapq_distribution.tsv"
out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1B.MAPQ_distribution/mapq_distribution.pdf"


mapq <- read.csv(file, col.names = c("sample", "mapq"), sep = "\t")

pC <- ggplot(mapq, aes(x = mapq, y = sample, fill = sample)) +
  geom_violin(scale = "width", trim = TRUE) +
  labs(
    x = "MAPQ",
    y = NULL,
    title = "Mapping quality distribution"
  ) +
  theme_classic() +
  theme(legend.position = "none")
pC

ggsave(out_file, pC)
