# This script calculates the reproducibility of methylation levels at the site level for 5mC in non-CpG context 
# and 6mA. It takes as input the output of script 06.01, which are TXT files with the following columns: 
# chr, start, end, methylation (proportion), coverage (depth). 
# It computes Pearson correlations between replicates at the site level.

library(data.table) # version 1.18.0
library(dplyr) # version 1.1.4


# ============================================================
# 1. Set parameters and locate files
# ============================================================

mod <- "CA"

home <- path.expand("~")
path_write <- file.path(home, "Reproducibility", "Site-level", mod, strand) # Path to write output files
dir.create(path_write, recursive = TRUE, showWarnings = FALSE)

dir <- file.path(home, "Data_methylation", "datasets_proportions", mod) # Output of script 06.02 and 06.03

files <- list.files(
  path = dir,
  pattern = paste0(mod, "_", strand),
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
      "methylation",
      "coverage",
      "strand"
    )
  )
  
  # Create unique genomic site identifier
  x[, site_id := paste(
    chr,
    start,
    end,
    strand,
    sep = ":"
  )]
  
  # Keep only variables needed for downstream analyses
  x <- x[, .(
    site_id,
    methylation,
    coverage
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
  function(x) anyDuplicated(x$site_id)
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
    function(x) x$site_id
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
    x[site_id %in% common_sites]
  }
)

# Check
sapply(
  methylation_common,
  nrow
)

length(common_sites)


# ============================================================
# 9. Keep site_id + methylation and rename methylation
#     according to replicate
# ============================================================

methylation_matrix <- lapply(
  seq_along(methylation_common),
  function(i) {
    
    x <- methylation_common[[i]][
      ,
      .(
        site_id,
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
      by = "site_id",
      all = FALSE
    )
    
  },
  methylation_matrix
)


# ============================================================
# 11. Calculate pairwise Pearson correlations
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

