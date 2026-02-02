
file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1D.post_filter/dataset.txt"

ds <- read.csv(file, header = T, sep = "\t")




# Exclude chromosome NC_002333.2
X_filtered <- ds %>%
  filter(X.rname != "NC_002333.2")

# Boxplot of mean depth per chromosome
ggplot(X_filtered, aes(x = X.rname, y = meandepth)) +
  geom_boxplot(fill = "steelblue", outlier.color = "red") +
  theme_minimal() +
  labs(x = "Chromosome", y = "Mean Depth",
       title = "Mean Sequencing Depth per Chromosome Across Replicates") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





