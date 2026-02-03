library(data.table)
library(dplyr)

regions_name <- "promoters"

# 1. Create dataset with whole-brain RRBS

samples_wb <- c("M1", "M2", "F1", "F2")
path_wb <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Chaterjee/methylation_regions/", regions_name)

ds_wb <- NULL   # initialize

for (sam in samples_wb) {
  
  # Read file
  file <- file.path(path_wb, paste0(sam, "_", regions_name, ".txt"))
  ds <- read.csv(file, sep = "\t", header = F)
  
  # Remove regions without calls
  ds <- ds[ds$V7 != 0, ]
  
  # Reduce dataset. 5th column is methylation percentage
  ds_red <- ds[, c(1:3,5)]
  
  # Percentage over 1
  ds_red$V5 <- as.numeric(ds_red$V5)
  ds_red$V5 <- ds_red$V5 / 100 
  
  # Include sample name in % column
  colnames(ds_red) <- c("chr", "start", "end", paste0(sam, "_wb_mpct"))
  
  # Merge into dataset
  if (is.null(ds_wb)) {
    # First iteration → initialize
    ds_wb <- ds_red
  } else {
    # Subsequent iterations → merge
    ds_wb <- merge(ds_wb, ds_red, by = c("chr", "start", "end"), all = FALSE)
  }
}


# 2. Create dataset with forebrain ONT 

samples_fb <- c("02", "03", "04", "05")
path_fb <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Chaterjee/forebrain_datasets/methylation_regions/", regions_name)

ds_fb <- NULL   # initialize

for (sam in samples_fb) {
  #sam <- "02"
  # Read file
  file <- file.path(path_fb, paste0(sam, "_", regions_name, ".txt"))
  ds <- read.csv(file, sep = "\t", header = F)
  
  # Remove regions without calls
  ds <- ds[ds$V7 != 0, ]
  
  # Reduce dataset. 5th column is methylation percentage
  ds_red <- ds[, c(1:3,5)]
  
  # Include sample name in % column
  colnames(ds_red) <- c("chr", "start", "end", paste0(sam, "_fb_mpct"))
  
  # Merge into dataset
  if (is.null(ds_fb)) {
    # First iteration → initialize
    ds_fb <- ds_red
  } else {
    # Subsequent iterations → merge
    ds_fb <- merge(ds_fb, ds_red, by = c("chr", "start", "end"), all = FALSE)
  }
}


# Join both large datasets
final_ds <- merge(ds_wb, ds_fb, by = c("chr", "start", "end"))




mat <- as.matrix(final_ds[, 4:ncol(final_ds)])
storage.mode(mat) <- "numeric"

pca <- prcomp(
  t(mat),
  center = TRUE,
  scale. = F
)

plot(pca$x[,1], pca$x[,2],
     xlab = "PC1",
     ylab = "PC2",
     pch = 19)


library(ggplot2)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  sample = rownames(pca$x)
)

ggplot(pca_df, aes(PC1, PC2, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.6) +
  theme_classic()




