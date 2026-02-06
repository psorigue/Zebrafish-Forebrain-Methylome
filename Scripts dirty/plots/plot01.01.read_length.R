library(ggplot2)
library(scales)

path <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1A.Read_length/"
out_file <- paste0(path, "read_length.pdf")

ds <- read.csv(paste0(path, "read_length_rep01.txt"), header = F, col.names = "read_length")


ds$replicate <- "rep1"

# Read files into one dataset
for (num in c("2", "3", "4", "5", "6")) {
  
  new_ds <- read.csv(paste0(path, "read_length_rep0", num,  ".txt"), header = F, col.names = "read_length")
  new_ds$replicate <- paste0("rep", num)
  ds <- rbind(ds, new_ds)
  
}
# Try with filtered
ds <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/QC/read_length_filtered/rep1.txt", header = F, col.names = "read_length")
ds$replicate <- "rep1"
for (num in c("2", "3", "4", "5", "6")) {
  
  new_ds <- read.csv(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/QC/read_length_filtered/rep", num,  ".txt"), header = F, col.names = "read_length")
  new_ds$replicate <- paste0("rep", num)
  ds <- rbind(ds, new_ds)
  
}

ds$replicate <- factor(ds$replicate, levels = paste0("rep", 6:1))

blue_reps <- c(
  "#dce9f5",
  "#bcd7ee",
  "#8fbde6",
  "#5fa2da",
  "#2f78c4",
  "#1f4e99"
)

p <- ggplot(ds, aes(x = read_length, fill = replicate)) +
  geom_histogram(bins = 100, position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = blue_reps) +
  scale_x_log10(
    labels = label_number(scale_cut = cut_si(unit = ""))
  ) +
  labs(
    x = "Read length (bp, log10)",
    y = "Read count",
    title = "Read length distribution",
    fill = "Replicate"
  ) +
  theme_classic()
p
b <- ggplot(ds, aes(x = read_length, y = replicate, fill = replicate)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_x_log10()
b
ggsave(out_file, b)
