path <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1A.Read_length/"


ds <- read.csv(paste0(path, "read_length_rep01.txt"), header = F, col.names = "read_length")
ds$replicate <- "rep1"

# Read files into one dataset
for (num in c("2", "3", "4", "5", "6")) {
  
  new_ds <- read.csv(paste0(path, "read_length_rep0", num,  ".txt"), header = F, col.names = "read_length")
  new_ds$replicate <- paste0("rep", num)
  ds <- rbind(ds, new_ds)
  
}



ggplot(ds, aes(x = read_length, fill = replicate)) +
  geom_histogram(bins = 100, position = "stack") +
  scale_x_log10(
    labels = label_number(scale_cut = cut_si(unit = ""))
  ) +
  labs(
    x = "Read length (bp, log10)",
    y = "Read count",
    title = "Read length distribution"
  ) +
  theme_classic()


ggplot(ds, aes(x = read_length, y = replicate, fill = replicate)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_x_log10()
