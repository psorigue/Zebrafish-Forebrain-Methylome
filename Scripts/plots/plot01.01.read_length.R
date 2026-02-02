file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1A.Read_length/"

read_lengths <- read_tsv("read_lengths.tsv")

pA <- ggplot(read_lengths, aes(x = read_length, y = replicate, fill = replicate)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_x_log10(labels = label_number_si()) +
  labs(
    x = "Read length (bp, log10)",
    y = NULL,
    title = "Read length distribution"
  ) +
  theme_classic() +
  theme(legend.position = "none")