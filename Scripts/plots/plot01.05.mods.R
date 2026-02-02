library(patchwork) 
library(ggplot2)

file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/1E.mods/dataset.txt"

ds <- read.csv(file, header = T, sep = "\t")


# Filter for base A
df_A <- ds %>%
  filter(base == "A" & code %in% c("can", "a")) %>%
  mutate(code = factor(code),
         replicate = factor(replicate, levels = unique(replicate)))

plot_A <- ggplot(df_A, aes(x = code, y = pass_frac, fill = replicate)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(title = "Base A", x = "Category", y = "Pass Fraction") +
  scale_fill_brewer(palette = "Set2")

# Filter for base C
df_C <- ds %>%
  filter(base == "C" & code %in% c("can", "m", "h")) %>%
  mutate(code = factor(code),
         replicate = factor(replicate, levels = unique(replicate)))

plot_C <- ggplot(df_C, aes(x = code, y = pass_frac, fill = replicate)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(title = "Base C", x = "Category", y = "Pass Fraction") +
  scale_fill_brewer(palette = "Set2")

# Combine plots side by side
plot_A + plot_C

  
  
  