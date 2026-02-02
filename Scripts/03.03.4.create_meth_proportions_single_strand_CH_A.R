# Script to create methylation proportions files from raw counts files for 5mC and 5hmC datasets. 
# As input, it takes the output of script 03.02.pileup.sh: TXT files with the following columns: chr, start, end, N (depth), X (methylated reads), strand. In CH case, it also includes a motif column to be able to separate contexts (CA, CC, CT).

library(data.table) # v.1.18

# ===============================
# 1. Paths
# ===============================
mod <- "ch"
path_in <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_by_mod/", mod)
path_out <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_proportions/", mod, "_ss/")
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

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
  
  # Split by pattern
  pat_list <- split(dt, by = "pattern", keep.by = TRUE)
  
  # Process each context
  for(pat in names(pat_list)) {
    dt_ctx <- pat_list[[pat]]
    
    dt_out <- dt_ctx[, .(chr, start, end, mpct, N, strand)]

    # Base output name
    base_name <- tools::file_path_sans_ext(basename(f))
    
    # Output files
    outfile <- file.path(path_out, paste0(base_name, "_", pat, "_mpct.txt"))
    
    # Write only if non-empty
    fwrite(dt_out, outfile, sep="\t", col.names=FALSE, quote=FALSE)
    
  }
  
  rm(dt, dt_ctx, pat_list)
  gc()
}
