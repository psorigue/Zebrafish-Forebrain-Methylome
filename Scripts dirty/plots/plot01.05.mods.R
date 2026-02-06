library(patchwork) 
library(ggplot2)
library(ggbreak)
library(dplyr)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1E.mods/dataset.txt"
out_pdf <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1E.mods/mods.pdf"
out_tiff <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1E.mods/mods.tiff"

ds <- read.csv(file, header = T, sep = "\t")

# Replace mod code
replacements <- c("can" = "canonical", "a" = "6mA", "m" = "5mC", "h" = "5hmC")
ds <- ds %>%
  mutate(code = recode(code, !!!replacements))

blue_reps <- c(
  "#dce9f5",
  "#bcd7ee",
  "#8fbde6",
  "#5fa2da",
  "#2f78c4",
  "#1f4e99"
)

# Filter for base A
df_A <- ds %>%
  filter(base == "A" & code %in% c("canonical", "6mA")) %>%
  mutate(code = factor(code),
         replicate = factor(replicate, levels = unique(replicate)))
df_A_frac_mod <- df_A[df_A$code == "6mA",]
df_A_frac_can <- df_A[df_A$code == "canonical",]
median_mod <- median(df_A_frac_mod$pass_frac)
median_can <- median(df_A_frac_can$pass_frac)
se_x   <- sd(df_A_frac$pass_frac) / sqrt(length(df_A_frac$pass_frac))
c(mean_x, se_x)

df_A$code <- factor(df_A$code, levels = c("canonical", "6mA"))
plot_A <- ggplot(df_A, aes(x = code, y = pass_frac, fill = replicate)) +
  geom_col(
    position = position_dodge(width = 0.6),
    width = 0.5,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = blue_reps) +
  labs(
    title = "Base A",
    x = "Category",
    y = "Fraction"
  ) +  # remove fill label
  theme_classic() +
  theme(legend.position = "none") +  # hide legend
  geom_text(
    data = df_A %>% group_by(code) %>% summarize(pass_frac = median(pass_frac)),
    aes(x = code, y = pass_frac, label = round(pass_frac, 3)),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 2
  )
# Filter for base C
df_C <- ds %>%
  filter(base == "C" & code %in% c("canonical", "5mC", "5hmC")) %>%
  mutate(code = factor(code),
         replicate = factor(replicate, levels = unique(replicate)))
df_C_can <- df_C[df_C$code == "canonical",]
df_C_m <- df_C[df_C$code == "5mC",]
df_C_h <- df_C[df_C$code == "5hmC",]
mean_can <- median(df_C_can$pass_frac)

se_x   <- sd(df_C_frac$pass_frac) / sqrt(length(df_C_frac$pass_frac))
c(mean_x, se_x)

df_C$code <- factor(df_C$code, levels = c("canonical", "5mC", "5hmC"))
plot_C <- ggplot(df_C, aes(x = code, y = pass_frac, fill = replicate)) +
  geom_col(
    position = position_dodge(width = 0.6),
    width = 0.5,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = blue_reps) +
  labs(
    title = "Base C",
    x = "Category",
    y = NULL,
    fill = "Replicate"
  ) +
  geom_text(
    data = df_C %>% group_by(code) %>% summarize(pass_frac = median(pass_frac)),
    aes(x = code, y = pass_frac, label = round(pass_frac, 3)),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 2) +
  theme_classic()
plot_C 
  
# Combine plots side by side
p <- plot_A + plot_C
p
ggsave(out_file, p)

# Make tiff
{
  tiff(
    out_tiff,
    width = 5,
    height = 4,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p)
  dev.off()
  
}
