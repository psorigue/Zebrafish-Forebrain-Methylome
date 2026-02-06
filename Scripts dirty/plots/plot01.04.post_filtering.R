library(ggplot2)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1D.post_filter/dataset.txt"
out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1D.post_filter/post_filter.pdf"

ds <- read.csv(file, header = T, sep = "\t")


# Exclude chromosome NC_002333.2
X_filtered <- ds %>%
  filter(X.rname != "NC_002333.2")

# Boxplot of mean depth per chromosome
p <- ggplot(X_filtered, aes(x = X.rname, y = meandepth)) +
  geom_boxplot(fill = "#8fbde6", width = 0.5) +
  theme_minimal() +
  labs(x = "Chromosome", y = "Mean Depth",
       title = "Mean Sequencing Depth per Chromosome Across Replicates") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(out_file, p)




