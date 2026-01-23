library(data.table)

mod <- "5hmC"

setwd(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_by_mod/", mod, "/"))
out_dir <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/meth_proportions/", mod)
dir.create(out_dir, showWarnings = F, recursive = T)


# Process all .txt files in current directory
files <- list.files(pattern = "\\.txt$")


for (f in files) {
  cat("Processing", f, "...\n")
  
  dt <- fread(f, header = FALSE, col.names = c("chr", "start", "end", "N", "X", "strand"))
  
  dt <- dt[, c("chr", "start", "end", "N", "X")]
  
  # Filter coverage = or > than 5
  dt <- dt[N > 4]
  
  # Calculate methylation ratio
  dt[, mpct := X / N]
  
  dt_final <- dt[, .(chr, start, end, N, mpct)]
  
  outfile <- file.path(
    out_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_mpct.txt")
  )
  
  fwrite(dt_final, outfile, sep = "\t", col.names = F, row.names = F, quote = F)
  
  rm(dt, dt_final)
  gc()
}
