library(data.table)
library(dplyr)

mod <- "6mA"

setwd(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_by_mod/", mod, "/"))
out_dir <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/meth_proportions/", mod)
dir.create(out_dir, showWarnings = F, recursive = T)


# Process all .txt files in current directory
files <- list.files(pattern = "\\.txt$")


for (f in files) {
  cat("Processing", f, "...\n")
  
  # Read file
  dt <- fread(
    f,
    header = FALSE,
    col.names = c("chr", "start", "end", "N", "X", "strand")
  )
  
  # Filter rows with coverage > 1
  dt <- dt[N > 1]
  
  # Calculate methylation ratio
  dt[, mpct := X / N]
  
  # Split by strand and keep only BED4 columns: chr, start, end, mpct
  dt_pos <- dt[strand == "+", .(chr, start, end, mpct)]
  dt_neg <- dt[strand == "-", .(chr, start, end, mpct)]
  
  # Base output name (without extension)
  base_name <- tools::file_path_sans_ext(basename(f))
  
  # Output file names
  outfile_pos <- file.path(out_dir, paste0(base_name, "_mpct_pos.txt"))
  outfile_neg <- file.path(out_dir, paste0(base_name, "_mpct_neg.txt"))
  
  # Write files (only if non-empty)
  fwrite(dt_pos, outfile_pos, sep = "\t")
  fwrite(dt_neg, outfile_neg, sep = "\t")
  
  # Clean up
  rm(dt, dt_pos, dt_neg)
  gc()
}  
