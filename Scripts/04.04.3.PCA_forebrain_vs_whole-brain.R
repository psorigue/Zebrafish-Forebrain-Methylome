library(data.table)
library(dplyr)

regions_name <- "genome_50kb_bins"

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
  
  # Reduce dataset. 4th column is methylation percentage
  ds_red <- ds[, 1:4]
  
  # Include sample name in % column
  colnames(ds_red) <- c("chr", "start", "end", paste0(sam, "_wb_pct"))
  
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
# Because we compare to RRBS, we must add up the signal that 5mC and 5hmC produced in nanopore. RRBS signal includes 5mC 
# and 5hmC modifications

samples_fb <- c("02", "03", "04", "05")
path_fb <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/", regions_name)

ds_fb <- NULL   # initialize

for (sam in samples_fb) {
  #sam <- "01"
  # Read files
  file_m <- file.path(path_fb, paste0("5mC/meth_mean_5mC_rep", sam, "_", regions_name, ".txt"))
  file_h <- file.path(path_fb, paste0("5hmC/meth_mean_5hmC_rep", sam, "_", regions_name, ".txt"))
  ds_m <- fread(file_m, sep = "\t", header = F)
  ds_h <- fread(file_h, sep = "\t", header = F)
  
  # Concatenate files. We have to sum methylation calls, keeping coverage as it is
  ds <- rbind(ds_m, ds_h)
  
  # Group same region and sum modification calls
  ds_grouped <- ds %>%
    group_by(V1, V2, V3) %>%
    summarise(
      V4 = sum(V4, na.rm = TRUE),
      V7 = max(V7, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    select(V1, V2, V3, V4, V7)
  
  # Remove regions without calls
  ds_grouped <- ds_grouped[ds_grouped$V7 != 0, ]
  
  # Reduce dataset. 4th column is methylation percentage
  ds_red <- ds_grouped[, 1:4]
  
  # Include sample name in % column
  colnames(ds_red) <- c("chr", "start", "end", paste0(sam, "_fb_pct"))
  
  # Merge into dataset
  if (is.null(ds_fb)) {
    # First iteration → initialize
    ds_fb <- ds_red
  } else {
    # Subsequent iterations → merge
    ds_fb <- merge(ds_fb, ds_red, by = c("chr", "start", "end"), all = FALSE)
  }
}