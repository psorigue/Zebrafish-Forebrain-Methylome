library(dplyr)

mods <- c("CT", "CA", "CC")
an <- "genome_50kb_bins"
N <- 100

setwd(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/output/", 
             an, "/ch/"))

num <- c("1", "2", "3", "4", "5", "6")
strand <- c("pos", "neg")

# Store results for all motifs
all_mods_list <- lapply(mods, function(mod) {
  
  # Build all file names for this motif
  files <- expand.grid(strand = strand, num = num) %>%
    mutate(file = paste0("meth_mean_", mod, "_", strand,
                         "_rep", num, "_", an, ".txt"))
  
  # Read files and add metadata columns
  df_list <- lapply(seq_len(nrow(files)), function(i) {
    read.csv(files$file[i], sep = "\t", header = FALSE) %>%
      mutate(
        replicate = paste0("rep", files$num[i]),
        strand = files$strand[i],
        motif = mod
      )
  })
  
  # Combine + filter
  bind_rows(df_list) %>%
    filter(V7 > N)   # adjust column if needed
})

# Combine all motifs into a single dataset
df_final <- bind_rows(all_mods_list)

# Write single output file
write.table(df_final,
            paste0("CH_", an, "_stranded.txt"),
            quote = FALSE,
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)