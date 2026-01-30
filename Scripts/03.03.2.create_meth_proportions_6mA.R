
# Script to create methylation proportions files from raw counts files for 5mC and 5hmC datasets. 
# As input, it takes the output of script 03.02.pileup.sh: TXT files with the following columns: chr, start, end, N (depth), X (methylated reads), strand. 

# Load required libraries
library(data.table)

# Set working directory and output directory
setwd(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_by_mod/6mA/"))
out_dir <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/meth_proportions/6mA")
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
  dt <- dt[N > 4]

  # Exclude mitochondrial chromosome
  dt <- dt[chr != "NC_002333.2"]
  
  # Calculate methylation ratio
  dt[, mpct := X / N]
  
  # Split by strand and keep only BED4 columns: chr, start, end, mpct, strand
  dt_pos <- dt[strand == "+", .(chr, start, end, mpct, N, strand)]
  dt_neg <- dt[strand == "-", .(chr, start, end, mpct, N, strand)]
  
  # Base output name (without extension)
  base_name <- tools::file_path_sans_ext(basename(f))
  
  # Output file names keeping strand information
  outfile_pos <- file.path(out_dir, paste0(base_name, "_mpct_pos.txt"))
  outfile_neg <- file.path(out_dir, paste0(base_name, "_mpct_neg.txt"))
  
  # Write files
  fwrite(dt_pos, outfile_pos, col.names = F, sep = "\t")
  fwrite(dt_neg, outfile_neg, col.names = F, sep = "\t")

  # Clean up memory
  rm(dt, dt_pos, dt_neg)
  gc()
}  
