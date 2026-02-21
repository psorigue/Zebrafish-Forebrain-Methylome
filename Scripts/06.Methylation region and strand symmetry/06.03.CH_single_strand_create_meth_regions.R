# Script to create methylation proportions files from raw counts files for 5mC modifications in CH motifs. 
# As input, it takes the output of script 02.pileup.sh: TXT files with the following columns: chr, start, end, N (depth), X (methylated reads), strand, and a motif column to be able to separate contexts (CA, CC, CT).
# This script creates single-strand files instead of separate files for each strand, to compute overall methylation levels in later analyses.

library(data.table) # v.1.18

# ===============================
# 1. Paths
# ===============================
home <- path.expand("~")
path_in <- paste0(home, "/Pol/Methylome/Data_methylation/datasets_by_mod/ch")
path_out <- paste0(home, "/Pol/Methylome/Data_methylation/datasets_proportions/ch_ss/")
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

# Define input files by full path
files <- list.files(path_in, pattern="\\.txt$", full.names=TRUE)

# ===============================
# 2. Process each file
# ===============================
for(f in files) {
  cat("Processing", basename(f), "...\n")
  
  # Read file (7 columns)
  dt <- fread(f,
              header = FALSE,
              col.names = c("chr", "start", "end", "N", "X", "strand", "motif"))
  
  # Extract dinucleotide pattern (CA, CC, CT) from motif column
  dt[, pattern := tstrsplit(motif, ",")[[2]]]
  
  # Filter coverage
  dt <- dt[N > 4]
  
  # Calculate methylation ratio
  dt[, mpct := X / N]
  
  # Exclude mitochondrial chromosome
  dt <- dt[chr != "NC_002333.2"]
  
  # Split by pattern column
  pat_list <- split(dt, by = "pattern", keep.by = TRUE)
  
  # Process each context
  for(pat in names(pat_list)) {
    dt_ctx <- pat_list[[pat]]
    
    # Combine strands into single file
    dt_out <- dt_ctx[, .(chr, start, end, mpct, N, strand)]

    # Base output name
    base_name <- tools::file_path_sans_ext(basename(f))
    
    # Output files
    outfile <- file.path(path_out, paste0(base_name, "_", pat, "_mpct.txt"))
    
    # Write only if non-empty
    fwrite(dt_out, outfile, sep="\t", col.names=FALSE, quote=FALSE)
    
  }
  
  # Free memory
  rm(dt, dt_ctx, pat_list)
  gc()
}
