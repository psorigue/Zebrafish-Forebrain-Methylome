library(dplyr)
library(ggplot2)

path_in <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/clustering/datasets/"
path_out <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/clustering/results/"

# Create output folder if it doesn't exist
dir.create(path_out, showWarnings = FALSE, recursive = TRUE)

for (motif in c("CA", "CC", "CT")) {
  
  ratio_list <- list()
  
  # Read each replicate and compute ratio
  for (rep in 1:6) {
    
    file_total <- paste0(path_in, motif, "_rep", rep, "_total.txt")
    file_filt  <- paste0(path_in, motif, "_rep", rep, "_filt.txt")
    
    ds_total <- read.csv(file_total, sep = "\t", header = F,
                         col.names = c("chr", "start", "end", "name", "count_total"))
    ds_filt  <- read.csv(file_filt, sep = "\t", header = F,
                         col.names = c("chr", "start", "end", "name", "count_filt"))
    
    ds_ratio <- merge(ds_total[, c("name", "count_total")],
                      ds_filt[, c("name", "count_filt")],
                      by = "name", all = F) %>%
      mutate(!!paste0("ratio_rep", rep) := count_filt / count_total) %>%
      select(name, starts_with("ratio_rep"))
    
    ratio_list[[rep]] <- ds_ratio
  }
  
  # Merge all replicates into one dataset
  all_ratios <- Reduce(function(x, y) merge(x, y, by = "name"), ratio_list)
  
  all_ratios <- na.omit(all_ratios)
  
  # Compute median across replicates
  all_ratios <- all_ratios %>%
    rowwise() %>%
    mutate(median_ratio = median(c_across(starts_with("ratio_rep")))) %>%
    ungroup()
  
  # Save the dataset with replicates + median
  write.csv(all_ratios, file = paste0(path_out, motif, "_ratios_median.csv"),
            row.names = FALSE)
  
  # Plot histogram of median ratios and save
  p <- ggplot(all_ratios, aes(x = median_ratio)) +
    geom_histogram(binwidth = 0.0005, fill = "skyblue", color = "black") +
    ggtitle(paste0("Median ratio per region for motif ", motif)) +
    xlab("Median ratio") +
    ylab("Number of regions") +
    theme_minimal()
  p
  ggsave(filename = paste0(path_out, motif, "_median_ratio_histogram.pdf"),
         plot = p, width = 6, height = 4)
}
