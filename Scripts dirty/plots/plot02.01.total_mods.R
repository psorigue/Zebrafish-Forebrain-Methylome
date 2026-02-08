library(dplyr)
file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/writing/ForebrainMethylome_paper/Data/3.modifications_motif_proportions/output_all.txt"

ds <- read.csv(file, header = T, sep = "\t")

ds_mod <- ds %>%
  mutate(
    mod_code = recode(
      mod_code,
      "m" = "5mC",
      "h" = "5hmC",
      "a" = "6mA"
    ))

# CpG
# Subset for CG motif and codes m and h
df_CG <- ds_mod %>%
  filter(motif == "CG", mod_code %in% c("5mC", "5hmC"))

# Calculate fraction mean and SE
df_CG_frac <- df_CG[df_CG$mod_code == "h",]
mean_x <- mean(df_CG_frac$frac_mod)
se_x   <- sd(df_CG_frac$frac_mod) / sqrt(length(df_CG_frac$frac_mod))
c(mean_x, se_x)

# Calculate median to include in the plot
df_CG_median <- df_CG %>%
  group_by(mod_code) %>%
  summarize(
    median_mod = median(mod_over_total),
    .groups = "drop"
  )

# Refactor
df_CG$replicate <- factor(df_CG$replicate, levels = unique(df_CG$replicate))
df_CG$mod_code <- factor(df_CG$mod_code, levels = c("5mC", "5hmC"))

# Bar plot
blue_reps <- c(
  "#dce9f5",
  "#bcd7ee",
  "#8fbde6",
  "#5fa2da",
  "#2f78c4",
  "#1f4e99"
)

p <- ggplot(df_CG, aes(x = mod_code, y = mod_over_total, fill = replicate)) +
  geom_col(
    position = position_dodge(width = 0.6),
    width = 0.5,
    color = "black",
    linewidth = 0.2
  ) +
  geom_text(
    data = df_CG_median,
    aes(
      x = mod_code,
      y = median_mod,
      label = round(median_mod, 4)
    ),
    inherit.aes = FALSE,
    vjust = -1,
    size = 3
  ) +
  scale_fill_manual(values = blue_reps) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "Modification",
    y = "Fraction Modified",
    title = "Fraction of Modification for CG Motifs Across Replicates"
  )
p

# nonCpG (only m)
# Subset for motifs CA, CC, CT and code m only
df_m <- ds_mod %>%
  filter(motif %in% c("CA", "CC", "CT"), mod_code == "5mC")

# Calculate fraction mean and SE
df_CH_frac <- df_m[df_m$motif == "CT",]
mean_x <- mean(df_CH_frac$mod_over_total)
se_x   <- sd(df_CH_frac$mod_over_total) / sqrt(length(df_CH_frac$mod_over_total))
c(mean_x, se_x)

# Make motif a factor for plotting order
df_m$motif <- factor(df_m$motif, levels = c("CA", "CT", "CC"))

# Calculate medians to include in plot:
df_m_median <- df_m %>%
  group_by(motif) %>%
  summarize(median_mod = round(median(mod_over_total), 4), .groups = "drop")

# Boxplot: one box per motif, all replicates included as points
b <- ggplot(df_m, aes(x = motif, y = mod_over_total, fill = replicate)) +
  geom_col(
    position = position_dodge(width = 0.6),
    width = 0.5,
    color = "black",
    linewidth = 0.2
  ) +
  geom_text(
    data = df_m_median,
    aes(x = motif, y = median_mod, label = median_mod),
    inherit.aes = FALSE,
    vjust = -3,
    size = 3
  ) +
  scale_fill_manual(values = blue_reps) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # headroom for labels
  theme_minimal() +
  labs(
    x = "Motif",
    y = "Fraction Modified",
    title = "Fraction of m Modification for CA, CC, CT Motifs Across Replicates"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
b
out_file2 <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/2A.total_mods/noncpg.pdf"
ggsave(out_file2, b)


library(patchwork)
p <- p + b
p

# A
df_a <- ds_mod[ds_mod$motif == "A",]
mean_x <- mean(df_a$mod_over_total)
se_x   <- sd(df_CH_frac$mod_over_total) / sqrt(length(df_CH_frac$mod_over_total))
c(mean_x, se_x)

# A
# Subset for motifA
df_a <- ds_mod %>%
  filter(motif %in% "A", mod_code == "6mA")

# Calculate medians to include in plot:
df_a_median <- df_a %>%
  group_by(motif) %>%
  summarize(median_mod = median(mod_over_total), .groups = "drop")

# Boxplot: one box per motif, all replicates included as points
a <- ggplot(df_a, aes(x = motif, y = mod_over_total, fill = replicate)) +
  geom_col(
    position = position_dodge(width = 0.6),
    width = 0.5,
    color = "black",
    linewidth = 0.2
  ) +
  geom_text(
    data = df_a_median,
    aes(x = motif, y = median_mod, label = median_mod),
    inherit.aes = FALSE,
    vjust = -0.5,  # smaller negative so it sits just above the bar
    size = 3
  ) +
  scale_fill_manual(values = blue_reps) +
  scale_y_continuous(
    limits = c(0, 0.003),           # set y-axis limits
    expand = expansion(mult = c(0, 0.05))  # small padding at top
  ) +
  theme_minimal() +
  labs(
    x = "Motif",
    y = "Fraction Modified",
    title = "Fraction of a Modification for A. Motifs Across Replicates"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  )
a
