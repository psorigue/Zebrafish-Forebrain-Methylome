

# Concatenate files adding replicate column
mod <- "CA"
an <- "genome_50kb_bins"
N <- 10


setwd(paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/output/", an, "/", mod, "_ss/"))

num <- c("1", "2", "3", "4", "5", "6")
files <- paste0("meth_mean_", mod, "_rep", num, "_", an, ".txt")

# Read and add replicate column
df_list <- lapply(seq_along(files), function(i) {
  read.csv(files[i], sep = "\t", header = F) %>% mutate(replicate = paste0("rep", i))
})

# Combine all
df_all <- bind_rows(df_list)

# Remove regions with less than N occurrences (sites covered)
df_fin <- df_all[df_all$V7 > N,] # col 7 or 9 depends if there is strandedness in regions



write.table(df_fin, paste0(mod, "_all_", an, ".txt"), quote = F, sep = "\t", col.names = F, row.names = F)
