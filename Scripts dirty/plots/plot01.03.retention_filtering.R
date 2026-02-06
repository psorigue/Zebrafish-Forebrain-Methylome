library(scales)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1C.filtering/dataset.txt"
out_pdf <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1C.filtering/filtering.pdf"
out_tiff <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1C.filtering/filtering.tiff"


filtering <- read.csv(file, header = TRUE, sep = "\t")
filtering$stage <- factor(filtering$stage, levels = c("initial", "filtered"))

pE <- ggplot(filtering, aes(x = replicate, y = read_count, fill = stage)) +
  geom_col(
    position = position_dodge(width = 0.6),
    width = 0.5,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "initial"  = "grey80",
      "filtered" = "#8fbde6"
    )
  ) +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_si(unit = ""))
  ) +
  labs(
    x = NULL,
    y = "Reads",
    fill = NULL,
    title = "Read retention after filtering"
  ) +
  theme_classic()


pE

ggsave(out_pdf, pE)

# Make tiff
{
  tiff(
    out_tiff,
    width = 7,
    height = 5,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(pE)
  dev.off()
  
}
