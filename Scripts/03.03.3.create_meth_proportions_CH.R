library(data.table)

# ===============================
# 1. Paths
# ===============================
mod <- "ch"
path_in <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Epi/Data_processed/datasets_by_mod/", mod)
path_out <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Epi/meth_proportions/", mod, "/")
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
  
  # Extract dinucleotide pattern (CA, CC, CT) from motif
  dt[, pattern := tstrsplit(motif, ",")[[2]]]
  
  # Filter coverage
  dt <- dt[N > 2]
  
  # Calculate methylation ratio
  dt[, mpct := X / N]
  
  # Split by pattern
  pat_list <- split(dt, by = "pattern", keep.by = TRUE)
  
  # Process each context
  for(pat in names(pat_list)) {
    dt_ctx <- pat_list[[pat]]
    
    # Split by strand
    dt_pos <- dt_ctx[strand == "+", .(chr, start, end, N, mpct)]
    dt_neg <- dt_ctx[strand == "-", .(chr, start, end, N, mpct)]
    
    # Base output name
    base_name <- tools::file_path_sans_ext(basename(f))
    
    # Output files
    outfile_pos <- file.path(path_out, paste0(base_name, "_", pat, "_pos.txt"))
    outfile_neg <- file.path(path_out, paste0(base_name, "_", pat, "_neg.txt"))
    
    # Write only if non-empty
    if(nrow(dt_pos) > 0) fwrite(dt_pos, outfile_pos, sep="\t", col.names=FALSE, quote=FALSE)
    if(nrow(dt_neg) > 0) fwrite(dt_neg, outfile_neg, sep="\t", col.names=FALSE, quote=FALSE)
    
    # Clean up
    rm(dt_pos, dt_neg)
  }
  
  rm(dt, dt_ctx, pat_list)
  gc()
}