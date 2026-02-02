

# Remove postive strand rows from "neg" file and opposite

# Reduce columns

# Plot

# For regions not needing strand: make this columns: chr start end region_id sample replicate mean_meth mean_cov motif

Minimum solid set:
  
  Density of mean methylation (per modification)

Violin/boxplot per motif

Replicate scatter

Coverage vs methylation

That already tells a coherent story.


# 1. Strand:
ggplot(df, aes(x = strand, y = mean_meth)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1) +
  theme_bw()

# 2. Coverage vs methylation
ggplot(df, aes(mean_cov, mean_meth)) +
  geom_point(alpha = 0.3) +
  theme_bw()

# 3. Coverage distribution
ggplot(df, aes(mean_cov)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  theme_bw()

# Methylation distribution
ggplot(df, aes(mean_meth)) +
  geom_density() +
  theme_bw()

# Replicate concordance: make columns