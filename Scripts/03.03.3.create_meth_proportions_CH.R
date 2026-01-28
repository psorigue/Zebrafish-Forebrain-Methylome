# Script to create methylation proportions files from raw counts files for 5mC and 5hmC datasets. 
# As input, it takes the output of script 03.02.pileup.sh: TXT files with the following columns: chr, start, end, N (depth), X (methylated reads), strand. In CH case, it also includes a motif column to be able to separate contexts (CA, CC, CT).

library(data.table) # v.1.18

# ===============================
# 1. Paths
# ===============================
mod <- "ch"
path_in <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_by_mod/", mod)
path_out <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/meth_proportions/", mod, "/")
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
  
  # Split by pattern
  pat_list <- split(dt, by = "pattern", keep.by = TRUE)
  
  # Process each context
  for(pat in names(pat_list)) {
    dt_ctx <- pat_list[[pat]]
    
    # Split by strand
    dt_pos <- dt_ctx[strand == "+", .(chr, start, end, mpct, N)]
    dt_neg <- dt_ctx[strand == "-", .(chr, start, end, mpct, N)]
    
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