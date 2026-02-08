library(scales)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1C.filtering/dataset.txt"
out_pdf <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1C.filtering/filtering.pdf"
out_tiff <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1C.filtering/filtering.tiff"


filtering <- read.csv(file, header = TRUE, sep = "\t")
filtering$stage <- factor(filtering$stage, levels = c("initial", "filtered"))

filtering$fill_group <- ifelse(
  filtering$stage == "filtered",
  as.character(filtering$replicate),
  "initial"
)

blue_reps <- c(
  "rep1" = "#dce9f5",
  "rep2" = "#bcd7ee",
  "rep3" = "#8fbde6",
  "rep4" = "#5fa2da",
  "rep5" = "#2f78c4",
  "rep6" = "#1f4e99"
)


pE <- ggplot(
  filtering,
  aes(x = replicate, y = read_count, fill = fill_group)
) +
  geom_col(
    position = position_dodge(width = 0.6),
    width = 0.5,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "initial" = "grey80",
      blue_reps   # named vector: replicate → blue shade
    )
  ) +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_si(unit = ""))
  ) +
  labs(
    x = "replicate",
    y = "Reads",
    fill = NULL,
    title = "Read retention after filtering"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14)
  )
pE


ggsave(out_pdf, pE)

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
  print(pE)
  dev.off()
  
}


