# This script reads the methylation data for each replicate, filters the bins by number of occurrences, and keeps only the common sites across all replicates. It then merges the data into a single matrix and calculates pairwise Pearson correlations between replicates at the region level for CpG sites.

library(data.table) # version
library(dplyr) # version 1.1.4


# ============================================================
# 1. Set parameters and locate files
# ============================================================

mod <- "5mC"

home <- path.expand("~")
path_write <- file.path(home, "Reproducibility", "Region-level", mod) # Path to write output files
dir.create(path_write, recursive = TRUE, showWarnings = FALSE)

dir <- file.path(home, "Data_methylation", "methylation_regions", "output", "genome_50kb_bins", mod)

# List all files in the directory with .txt extension
files <- list.files(
  path = dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

# Check files
files
basename(files)


# ============================================================
# 2. Define clean replicate names
# ============================================================

# Use the file names as replicate names
replicate_numbers <- c("1", "2", "3", "4", "5", "6")

replicate_names <- paste0("rep0", replicate_numbers)



# ============================================================
# 3. Function to read each methylation dataset
# ============================================================

read_methylation <- function(file) {
  
  x <- fread(
    file,
    header = FALSE,
    col.names = c(
      "chr",
      "start",
      "end",
      "binID",
      "methylation",
      "coverage",
      "num_motifs"
    )
  )
  
  # Filter bins by number of occurrences
  x <- x[
      num_motifs >= 100
  ]
  
  # Keep only variables needed for downstream analyses
  x <- x[, .(
    binID,
    methylation,
    coverage,
    num_motifs
  )]
  
  return(x)
}


# ============================================================
# 4. Load all six datasets
# ============================================================

methylation <- lapply(
  files,
  read_methylation
)

names(methylation) <- replicate_names


# ============================================================
# 5. Check for duplicated sites
# ============================================================

sapply(
  methylation,
  function(x) anyDuplicated(x$binID)
)


# ============================================================
# 6. Check number of sites per replicate
# ============================================================

# Total number of sites
a <- sapply(
  methylation,
  nrow
)

write.table(a, file = file.path(path_write, "sites_per_replicate.txt"), sep = "\t", quote = F)

# ============================================================
# 7. Find sites shared by ALL six replicates
# ============================================================

common_sites <- Reduce(
  intersect,
  lapply(
    methylation,
    function(x) x$binID
  )
)

b <- length(common_sites)

write.table(b, file = file.path(path_write, "number_common_sites.txt"), sep = "\t", quote = F)

# Proportion of each replicate retained
c <- sapply(
  methylation,
  function(x) {
    length(common_sites) / nrow(x)
  }
)

write.table(c, file = file.path(path_write, "proportion_retained_replicate.txt"), sep = "\t", quote = F)

# ============================================================
# 8. Keep only common sites
# ============================================================

methylation_common <- lapply(
  methylation,
  function(x) {
    x[binID %in% common_sites]
  }
)

# Check
sapply(
  methylation_common,
  nrow
)

length(common_sites)


# ============================================================
# 9. Keep binID + methylation and rename methylation
#     according to replicate
# ============================================================

methylation_matrix <- lapply(
  seq_along(methylation_common),
  function(i) {
    
    x <- methylation_common[[i]][
      ,
      .(
        binID,
        methylation
      )
    ]
    
    # Rename methylation column to replicate name
    setnames(
      x,
      "methylation",
      names(methylation_common)[i]
    )
    
    return(x)
  }
)



# ============================================================
# 10. Merge the six replicates
# ============================================================

methylation_matrix <- Reduce(
  function(x, y) {
    
    merge(
      x,
      y,
      by = "binID",
      all = FALSE
    )
    
  },
  methylation_matrix
)

# make sure methylation values are numeric
methylation_matrix[, (names(methylation_matrix)[-1]) := 
                     lapply(.SD, as.numeric),
                   .SDcols = names(methylation_matrix)[-1]]



# ============================================================
# 11. Convert to long format for plotting
# ============================================================

methylation_long <- melt(
  methylation_matrix,
  id.vars = "binID",
  variable.name = "replicate",
  value.name = "methylation"
)


# ============================================================
# 12. Calculate pairwise Pearson correlations
# ============================================================

methylation_values <- as.matrix(
  methylation_matrix[, -1]
)

pearson_cor <- cor(
  methylation_values,
  method = "pearson",
  use = "pairwise.complete.obs"
)

write.table(pearson_cor, file = file.path(path_write, "pearson_matrix.txt"), sep = "\t", quote = F)
