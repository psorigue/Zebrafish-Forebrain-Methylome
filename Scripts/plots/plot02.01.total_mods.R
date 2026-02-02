
file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2A.total_mods//dataset.txt"

ds <- read.csv(file, header = T, sep = "\t")


# CpG
# Subset for CG motif and codes m and h
df_CG <- ds %>%
  filter(motif == "CG", mod_code %in% c("m", "h"))

# Optional: make 'replicate' a factor for proper ordering
df_CG$replicate <- factor(df_CG$replicate, levels = unique(df_CG$replicate))

# Boxplot
ggplot(df_CG, aes(x = mod_code, y = frac_mod)) +
  geom_boxplot(fill = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = "Modification", y = "Fraction Modified", 
       title = "Fraction of Modification for CG Motifs Across Replicates") +
  scale_color_brewer(palette = "Set2")



# nonCpG (only m)
# Subset for motifs CA, CC, CT and code m only
df_m <- ds %>%
  filter(motif %in% c("CA", "CC", "CT"), mod_code == "m")

# Sum frac_mod per motif and replicate (usually one row per motif/replicate for m)
df_sum <- df_m %>%
  group_by(motif, replicate) %>%
  summarise(sum_frac = sum(frac_mod), .groups = "drop")

# Make motif a factor for plotting order
df_sum$motif <- factor(df_sum$motif, levels = c("CA", "CT", "CC"))

# Boxplot: one box per motif, all replicates included as points
ggplot(df_sum, aes(x = motif, y = sum_frac)) +
  geom_boxplot(fill = c("steelblue", "tomato", "gold")) +
  theme_minimal() +
  labs(x = "Motif", y = "Fraction Modified (m only)",
       title = "Fraction of m Modification for CA, CC, CT Motifs Across Replicates") +
  scale_color_brewer(palette = "Set2")
