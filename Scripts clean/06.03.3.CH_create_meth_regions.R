# Script to create methylation proportions files from raw counts files for 5mC modifications in CH motifs. 
# As input, it takes the output of script 02.pileup.sh: TXT files with the following columns: chr, start, end, N (depth), X (methylated reads), strand, and a motif column to be able to separate contexts (CA, CC, CT).
# Note that here we create separate files for each strand.

library(data.table) # v.1.18

# ===============================
# 1. Paths
# ===============================
path_in <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_by_mod/ch"
path_out <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_proportions/ch/"
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
  
  # Split by pattern
  pat_list <- split(dt, by = "pattern", keep.by = TRUE)
  
  # Process each context
  for(pat in names(pat_list)) {
    dt_ctx <- pat_list[[pat]]
    
    # Split by strand
    dt_pos <- dt_ctx[strand == "+", .(chr, start, end, mpct, N, strand)]
    dt_neg <- dt_ctx[strand == "-", .(chr, start, end, mpct, N, strand)]
    
    # Base output name
    base_name <- tools::file_path_sans_ext(basename(f))
    
    # Output files
    outfile_pos <- file.path(path_out, paste0(base_name, "_", pat, "_pos_mpct.txt"))
    outfile_neg <- file.path(path_out, paste0(base_name, "_", pat, "_neg_mpct.txt"))
    
    # Write only if non-empty
    if(nrow(dt_pos) > 0) fwrite(dt_pos, outfile_pos, sep="\t", col.names=FALSE, quote=FALSE)
    if(nrow(dt_neg) > 0) fwrite(dt_neg, outfile_neg, sep="\t", col.names=FALSE, quote=FALSE)
    
    # Clean up
    rm(dt_pos, dt_neg)
  }
  
  # Free memory
  rm(dt, dt_ctx, pat_list)
  gc()
}