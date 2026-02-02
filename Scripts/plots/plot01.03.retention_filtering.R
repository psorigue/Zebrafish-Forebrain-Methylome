library(scales)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1C.filtering/dataset.txt"

filtering <- read.csv(file, header = T, sep = "\t")

pE <- ggplot(filtering, aes(x = replicate, y = read_count, fill = stage)) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_si(unit = ""))  # <-- note unit=""
  ) +
  labs(
    x = NULL,
    y = "Reads",
    fill = NULL,
    title = "Read retention after filtering"
  ) +
  theme_classic()
pE
