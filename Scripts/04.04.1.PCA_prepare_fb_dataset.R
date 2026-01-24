

# This script joins sums the calls of 5mC and 5hmC of forebrain (ONT) dataset
# This is done because RRBS data outputs both modification in same calls

# We add the calls and create meth_proportions. Then following script is bedtools map


library(data.table)
library(dplyr)


path <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_by_mod/"
out_dir <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Chaterjee/forebrain_datasets/meth_proportions_cpg")
dir.create(out_dir, showWarnings = F, recursive = T)


samples <- c("02", "03", "04", "05")

for (sam in samples) {
  #sam <- "02"
  
  file_m <- paste0(path, "/5mC/5mC_sample", sam, ".txt")
  file_h <- paste0(path, "/5hmC/5hmC_sample", sam, ".txt")
  
  dt_m <- fread(file_m, header = FALSE, col.names = c("chr", "start", "end", "N", "X", "strand"))
  dt_h <- fread(file_h, header = FALSE, col.names = c("chr", "start", "end", "N", "X", "strand"))
  
  # Concatenate files. We have to sum methylation calls
  ds <- rbind(dt_m, dt_h)
  
  # Group same region and sum modification calls
  ds_grouped <- ds %>%
    group_by(chr, start, end) %>%
    summarise(
      N = max(N, na.rm = TRUE),
      X = sum(X, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(N > 4, chr != "NC_002333.2") %>% # Filter coverage = or > than 5 and Remove mitochondrial chromosome
    mutate(mpct = X / N) %>% # Calculate methylation ratio
    select(chr, start, end, N, mpct)
    
  
  outfile <- file.path(
    out_dir, paste0(sam, "_cpg.bed"))
  
  fwrite(ds_grouped, outfile, sep = "\t", col.names = F, row.names = F, quote = F)
  
  rm(ds, ds_grouped)
  gc()

}