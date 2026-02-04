
file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/writing/ForebrainMethylome_paper/Data/3.modifications_motif_proportions/output_all.txt"

ds <- read.csv(file, header = T, sep = "\t")

# CpG
# Subset for CG motif and codes m and h
df_CG <- ds %>%
  filter(motif == "CG", mod_code %in% c("m", "h"))

# Calculate fraction mean and SE
df_CG_frac <- df_CG[df_CG$mod_code == "h",]
mean_x <- mean(df_CG_frac$frac_mod)
se_x   <- sd(df_CG_frac$frac_mod) / sqrt(length(df_CG_frac$frac_mod))
c(mean_x, se_x)

# Optional: make 'replicate' a factor for proper ordering
df_CG$replicate <- factor(df_CG$replicate, levels = unique(df_CG$replicate))

# Boxplot
p <- ggplot(df_CG, aes(x = mod_code, y = mod_over_total)) +
  geom_boxplot(fill = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = "Modification", y = "Fraction Modified", 
       title = "Fraction of Modification for CG Motifs Across Replicates") +
  scale_color_brewer(palette = "Set2")
p
out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2A.total_mods/cpg.pdf"
ggsave(out_file, p)

# nonCpG (only m)
# Subset for motifs CA, CC, CT and code m only
df_m <- ds %>%
  filter(motif %in% c("CA", "CC", "CT"), mod_code == "m")


# Calculate fraction mean and SE
df_CH_frac <- df_m[df_m$motif == "CT",]
mean_x <- mean(df_CH_frac$mod_over_total)
se_x   <- sd(df_CH_frac$mod_over_total) / sqrt(length(df_CH_frac$mod_over_total))
c(mean_x, se_x)

# Sum frac_mod per motif and replicate (usually one row per motif/replicate for m)
df_sum <- df_m %>%
  group_by(motif, replicate) %>%
  summarise(sum_frac = sum(mod_over_total), .groups = "drop")

# Make motif a factor for plotting order
df_sum$motif <- factor(df_sum$motif, levels = c("CA", "CT", "CC"))

# Boxplot: one box per motif, all replicates included as points
b <- ggplot(df_sum, aes(x = motif, y = sum_frac)) +
  geom_boxplot(fill = c("steelblue", "tomato", "gold")) +
  theme_minimal() +
  labs(x = "Motif", y = "Fraction Modified (m only)",
       title = "Fraction of m Modification for CA, CC, CT Motifs Across Replicates") +
  scale_color_brewer(palette = "Set2")
b
out_file2 <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2A.total_mods/noncpg.pdf"
ggsave(out_file2, b)


# A
df_a <- ds[ds$motif == "A",]
mean_x <- mean(df_a$mod_over_total)
se_x   <- sd(df_CH_frac$mod_over_total) / sqrt(length(df_CH_frac$mod_over_total))
c(mean_x, se_x)
