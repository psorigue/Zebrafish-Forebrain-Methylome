
regions_name <- "genome_50kb_bins"

# 1. Create dataset with whole-brain RRBS

samples_wb <- c("M1", "M2", "F1", "F2")
path_wb <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Chaterjee/methylation_regions/", regions_name)

ds_wb <- NULL   # initialize

for (sam in samples) {
  
  # Read file
  file <- file.path(path_wb, paste0(sam, "_", regions_name, ".txt"))
  ds <- read.csv(file, sep = "\t", header = F)
  
  # Remove positions without coverage
  ds <- ds[ds$V6 != 0, ]
  
  # Reduce dataset. 4th column is methylation percentage
  ds_red <- ds[, 1:4]
  
  # Include sample name in % column
  colnames(ds_red) <- c("chr", "start", "end", paste0(sam, "_m_pct"))
  
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
path_fb <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_by_mod/", regions_name)

ds_fb <- NULL   # initialize

for (sam in samples) {
  
  # Read files
  file_m <- file.path(path_wb, paste0(sam, "_genome_50kb_bins.txt"))
  ds <- read.csv(file, sep = "\t", header = F)
  
  # Remove positions without coverage
  ds <- ds[ds$V6 != 0, ]
  
  # Reduce dataset. 4th column is methylation percentage
  ds_red <- ds[, 1:4]
  
  # Include sample name in % column
  colnames(ds_red) <- c("chr", "start", "end", paste0(sam, "_m_pct"))
  
  # Merge into dataset
  if (is.null(ds_wb)) {
    # First iteration → initialize
    ds_wb <- ds_red
  } else {
    # Subsequent iterations → merge
    ds_wb <- merge(ds_wb, ds_red, by = c("chr", "start", "end"), all = FALSE)
  }
}