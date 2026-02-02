
# Script to create methylation proportions files from raw counts files for 5mC and 5hmC datasets. 
# As input, it takes the output of script 03.02.pileup.sh: TXT files with the following columns: chr, start, end, N (depth), X (methylated reads), strand. 

# Load required library
library(data.table)

# Define modification type
mod <- "5hmC" # 5mC and 5hmC

# Set working directory and output directory
setwd(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_by_mod/", mod, "/"))
out_dir <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_proportions/", mod)
dir.create(out_dir, showWarnings = F, recursive = T)

# Process all .txt files in current directory
files <- list.files(pattern = "\\.txt$")

for (f in files) {
  cat("Processing", f, "...\n")
  
  # Read data
  dt <- fread(f, header = FALSE, col.names = c("chr", "start", "end", "N", "X", "strand"))
  
  # Keep only relevant columns. We discard strand, as we only have '+' in modkit output for CpG methylation.
  dt <- dt[, c("chr", "start", "end", "N", "X")]
  
  # Filter coverage
  dt <- dt[N > 4]

  # Exclude mitochondrial chromosome
  dt <- dt[chr != "NC_002333.2"]
  
  # Calculate methylation ratio
  dt[, mpct := X / N]
  
  # Keep only relevant columns for BED4 compatibility
  dt_final <- dt[, .(chr, start, end, mpct, N)]
  
  # Define output file name
  outfile <- file.path(
    out_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_mpct.txt")
  )
  # Write output
  fwrite(dt_final, outfile, col.names = F, sep = "\t")
  
  # Clean up memory
  rm(dt, dt_final)
  gc()
}