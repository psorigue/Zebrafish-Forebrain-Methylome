library(dplyr)
library(tidyr)


regions_name<- "cgi"
path_fb <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Chaterjee/methylation_cgi/forebrain/"
path_wb <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Chaterjee/methylation_cgi/whole-brain/"

# Threshold number sites per region (N) and mean coverage (C)
N <- 20
C <- 10

# 1. Load and process forebrain and whole-brain datasets. The output has columns: chr, start, end, name, mean_meth, mean_cov, nCpG
forebrain <- list.files(path_fb, full.names=TRUE) %>%
  lapply(function(f){
    df <- read.table(f, header=FALSE, col.names=c("chr","start","end", "name", "meth","cov","nCpG"))
    # Convert numeric columns explicitly
    df$start <- as.numeric(df$start)
    df$end <- as.numeric(df$end)
    df$meth <- as.numeric(df$meth) 
    df$cov <- as.numeric(df$cov)
    df$nCpG <- as.numeric(df$nCpG)
    df
  }) %>%
  bind_rows(.id="sample") %>%
  # keep only rows where meth and cov exist
  filter(!is.na(meth) & !is.na(cov)) %>%
  group_by(chr,start,end) %>%
  summarise(
    mean_meth = mean(meth, na.rm = TRUE),
    mean_cov  = mean(cov, na.rm = TRUE),
    nCpG      = sum(nCpG, na.rm = TRUE),
    name      = first(name),
    .groups   = "drop"
  )

wholebrain <- list.files(path_wb, full.names=TRUE) %>%
  lapply(function(f){
    df <- read.table(f, header=FALSE, col.names=c("chr","start","end", "name","meth","cov","nCpG"))
    # Convert numeric columns explicitly
    df$start <- as.numeric(df$start)
    df$end <- as.numeric(df$end)
    df$meth <- as.numeric(df$meth) / 100 # To match forebrain dataset (percentage to proportion)
    df$cov <- as.numeric(df$cov)
    df$nCpG <- as.numeric(df$nCpG)
    df
  }) %>%
  bind_rows(.id="sample") %>%
  # keep only rows where meth and cov exist
  filter(!is.na(meth) & !is.na(cov)) %>%
  group_by(chr,start,end) %>%
  summarise(
    mean_meth = mean(meth, na.rm=TRUE),
    mean_cov = mean(cov, na.rm=TRUE),
    nCpG = sum(nCpG, na.rm=TRUE),
    .groups="drop")


# 2. Merge datasets on the same regions
merged <- inner_join(forebrain, wholebrain,
                     by=c("chr","start","end"),
                     suffix=c("_FB","_WB"))

# 3. Keep only regions with at least N CpGs in both datasets and coverage threshold
merged_filtered <- merged %>%
  filter(nCpG_FB >= N & nCpG_WB >= N) %>%
  filter(mean_cov_FB >= C & mean_cov_WB >= C)

# OPTIONAL: See how threshold of number of sites affects
thresholds <- 1:30
cors <- sapply(thresholds, function(t){
  df <- merged %>% filter(nCpG_FB >= t & nCpG_WB >= t)
  cor(df$mean_meth_FB, df$mean_meth_WB, method="pearson")
})
plot(thresholds, cors, type="b", xlab="Min CpG per region", ylab="Pearson correlation")

# 4. Simple Pearson and Spearman
cor_pearson <- cor(merged_filtered$mean_meth_FB, merged_filtered$mean_meth_WB, method="pearson")
cor_spearman <- cor(merged_filtered$mean_meth_FB, merged_filtered$mean_meth_WB, method="spearman")

# 5. Scatter plot to visualize correlation
library(ggplot2)
p <- ggplot(merged_filtered, aes(x=mean_meth_WB, y=mean_meth_FB)) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", color="#1f4e99") +
  xlab("Whole-brain methylation") + ylab("Forebrain methylation")
p
out_tiff <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/3A.FB_vs_WB/linear_cgi_fb_vs_wb.tiff"
# Make tiff
{
  tiff(
    out_tiff,
    width = 6,
    height = 5,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p)
  dev.off()
  
}

# 6. Analyze directional differences in methylation
merged_filtered <- merged_filtered %>%
  mutate(delta = mean_meth_FB - mean_meth_WB)

# Classify changes
threshold <- 0.15 # Percentage of delta
merged_filtered <- merged_filtered %>%
  mutate(direction = case_when(
    delta > threshold ~ "Hyper in FB",
    delta < -threshold ~ "Hypo in FB",
    TRUE ~ "No change"
  ))
table(merged_filtered$direction)

# Histogram of delta methylation
merged_filtered$direction <- gsub(merged_filtered$direction, pattern = "Hypo in FB", replacement = "Hypomethylated in forebrain")
merged_filtered$direction <- gsub(merged_filtered$direction, pattern = "Hyper in FB", replacement = "Hypermethylated in forebrain")

p2 <- ggplot(merged_filtered, aes(x=delta, fill=direction)) +
  geom_histogram(binwidth=0.02, position="stack", color="black", alpha=0.8) +
  geom_vline(xintercept=0, linetype="dashed", color="darkgray") +
  scale_fill_manual(values=c("Hypomethylated in forebrain"="#1f4e99", "Hypermethylated in forebrain"="gold", "No change"="gray"),
                    labels = c(
                      "Hypomethylated\nin forebrain",
                      "Hypermethylated\nin forebrain",
                      "No change"
                    )) +
  theme_minimal() +
  labs(
    x = expression(Delta~Methylation),
    y = "Number of CG islands",
    fill = "Direction",
    title = NULL
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face="bold")
  )

p2
out_tiff2 <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/3A.FB_vs_WB/delta_cgi_fb_vs_wb.tiff"
# Make tiff
{
  tiff(
    out_tiff2,
    width = 6,
    height = 5,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p2)
  dev.off()
  
}

library(patchwork)
a <- p + p2
a
out_tiff3 <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/plot/3A.FB_vs_WB/panel.tiff"
# Make tiff
{
  tiff(
    out_tiff3,
    width = 12,
    height = 5,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(a)
  dev.off()
  
}




merged_filtered[order(abs(merged_filtered$delta), decreasing = T),]


# Order and write merged_filtered table to file
out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Chaterjee/cgi_FB_vs_WB.txt"

merged_filtered_ord <- merged_filtered[order(abs(merged_filtered$delta), decreasing = T),]
write.table(merged_filtered_ord, out_file, sep = "\t", quote = F, col.names = T, row.names = F)

