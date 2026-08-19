library(data.table)
library(dplyr)
library(ggplot2)

# ============================================================
# 1. Set parameters and locate files
# ============================================================

mod <- "6mA"
strand <- "pos"

path_write <- file.path("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/reproducibility/site-level/", mod, strand)
dir.create(path_write, recursive = T, showWarnings = F)

dir <- file.path(
  "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Data_methylation/datasets_proportions/", mod
)

files <- list.files(
  path = dir,
  pattern = paste0("_", strand),
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
# 4. Function to read each methylation dataset
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
# 5. Load all six datasets
# ============================================================

methylation <- lapply(
  files,
  read_methylation
)

names(methylation) <- replicate_names


# ============================================================
# 6. Check for duplicated sites
# ============================================================

sapply(
  methylation,
  function(x) anyDuplicated(x$site_id)
)


# ============================================================
# 7. Check number of sites per replicate
# ============================================================

# Total number of sites
a <- sapply(
  methylation,
  nrow
)

write.table(a, file = file.path(path_write, "sites_per_replicate.txt"), sep = "\t", quote = F)


# ============================================================
# 8. Find sites shared by ALL six replicates
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
# 9. Keep only common sites
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
# 10. Keep site_id + methylation and rename methylation
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
# 11. Merge the six replicates
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
# 12. Check the resulting matrix
# ============================================================

head(methylation_matrix)

dim(methylation_matrix)

colnames(methylation_matrix)


# Check for missing values
colSums(
  is.na(methylation_matrix)
)


# Check methylation values
summary(
  as.matrix(
    methylation_matrix[, -1]
  )
)


# ============================================================
# 13. Convert to long format for plotting
# ============================================================

methylation_long <- melt(
  methylation_matrix,
  id.vars = "site_id",
  variable.name = "replicate",
  value.name = "methylation"
)

head(methylation_long)


# ============================================================
# 14. Plot methylation distributions
# ============================================================

ggplot(
  methylation_long,
  aes(
    x = methylation,
    fill = replicate
  )
) +
  geom_density(
    alpha = 0.3
  ) +
  labs(
    x = "Methylation (%)",
    y = "Density",
    fill = "Replicate"
  ) +
  theme_classic()


# ============================================================
# 15. Calculate pairwise Pearson correlations
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



# ============================================================
# 17. Pearson correlation heatmap
# ============================================================

# Convert matrix to data frame first
pearson_df <- as.data.frame(pearson_cor)
pearson_df$Replicate1 <- rownames(pearson_df)

# Melt using data.table (or standard data.frame method)
pearson_long <- melt(
  as.data.table(pearson_df), 
  id.vars = "Replicate1", 
  variable.name = "Replicate2", 
  value.name = "Correlation"
)

ggplot(
  pearson_long,
  aes(
    x = Replicate1,
    y = Replicate2,
    fill = Correlation
  )
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.2f", Correlation)),
    size = 4
  ) +
  scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0,
    name = "Pearson r"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Site-level methylation reproducibility"
  ) +
  coord_equal() +
  theme_classic()


# ============================================================
# 19. Prepare data for PCA
# ============================================================

pca_data <- t(methylation_values)

dim(pca_data)

# ============================================================
# 20. Remove constant CpG sites
# ============================================================

# Calculate variance of each CpG across the six replicates
site_variance <- apply(
  pca_data,
  2,
  var
)

# Keep only CpGs with non-zero variance
pca_data_variable <- pca_data[
  ,
  site_variance > 0,
  drop = FALSE
]

# Check how many sites remain
dim(pca_data_variable)

# Number of constant sites removed
sum(site_variance == 0)

# ============================================================
# 20. Run PCA
# ============================================================

pca <- prcomp(
  pca_data_variable,
  center = TRUE,
  scale. = TRUE
)


# ============================================================
# 21. PCA variance explained
# ============================================================

pca_variance <- pca$sdev^2 /
  sum(pca$sdev^2)

pca_variance

round(
  pca_variance * 100,
  2
)

# ============================================================
# 22. PCA plot
# ============================================================

pca_scores <- as.data.frame(
  pca$x
)

pca_scores$replicate <- rownames(pca_scores)

ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    label = replicate
  )
) +
  geom_point(
    size = 4
  ) +
  geom_text(
    vjust = -0.8,
    size = 4
  ) +
  labs(
    x = paste0(
      "PC1 (",
      round(pca_variance[1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance[2] * 100, 1),
      "%)"
    ),
    title = "PCA of site-level methylation"
  ) +
  theme_classic()
