library(dplyr)
library(tidyr)


regions_name<- "genome_50kb_bins"
path_fb <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Chaterjee/forebrain_datasets/methylation_regions/", regions_name, "/")
path_wb <- paste0("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/Chaterjee/methylation_regions/", regions_name, "/")


# Example: combine all forebrain replicates
forebrain <- list.files(path_fb, full.names=TRUE) %>%
  lapply(function(f){
    df <- read.table(f, header=FALSE, col.names=c("chr","start","end", "name", "meth","cov","nCpG"))
    # Convert numeric columns explicitly
    df$start <- as.numeric(df$start)
    df$end <- as.numeric(df$end)
    df$meth <- as.numeric(df$meth) * 100 # To match wholebrain
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

wholebrain <- list.files(path_wb, full.names=TRUE) %>%
  lapply(function(f){
    df <- read.table(f, header=FALSE, col.names=c("chr","start","end", "name","meth","cov","nCpG"))
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
    mean_meth = mean(meth, na.rm=TRUE),
    mean_cov = mean(cov, na.rm=TRUE),
    nCpG = sum(nCpG, na.rm=TRUE),
    .groups="drop")



# Merge datasets on the same regions
merged <- inner_join(forebrain, wholebrain,
                     by=c("chr","start","end"),
                     suffix=c("_FB","_WB"))

# Keep only regions with at least N CpGs in both datasets and coverage threshold
N <- 100
C <- 10
merged_filtered <- merged %>%
  filter(nCpG_FB >= N & nCpG_WB >= N) %>%
  filter(mean_cov_FB > C & mean_cov_WB > C)

# See how threshold affects
thresholds <- 1:30
cors <- sapply(thresholds, function(t){
  df <- merged %>% filter(nCpG_FB >= t & nCpG_WB >= t)
  cor(df$mean_meth_FB, df$mean_meth_WB, method="pearson")
})
plot(thresholds, cors, type="b", xlab="Min CpG per region", ylab="Pearson correlation")

# Simple Pearson and Spearman
cor_pearson <- cor(merged_filtered$mean_meth_FB, merged_filtered$mean_meth_WB, method="pearson")
cor_spearman <- cor(merged_filtered$mean_meth_FB, merged_filtered$mean_meth_WB, method="spearman")


# Scatter plot
library(ggplot2)
ggplot(merged, aes(x=mean_meth_WB, y=mean_meth_FB)) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", color="red") +
  xlab("Whole-brain RRBS") + ylab("Forebrain ONT") +
  ggtitle(paste0("Pearson: ", round(cor_pearson,2),
                 " Spearman: ", round(cor_spearman,2)))

#install.packages("weights")
library(weights)
wtd.cor(merged$mean_meth_FB, merged$mean_meth_WB, weight=merged$nCpG_FB + merged$nCpG_WB)





# Directional differences
merged_filtered <- merged_filtered %>%
  mutate(delta = mean_meth_FB - mean_meth_WB)

# Classify changes
threshold <- 20 # Percentage of delta
merged_filtered <- merged_filtered %>%
  mutate(direction = case_when(
    delta > threshold ~ "Hyper in FB",
    delta < -threshold ~ "Hypo in FB",
    TRUE ~ "No change"
  ))
table(merged_filtered$direction)

# Histogram of delta methylation
ggplot(merged_filtered, aes(x=delta, fill=direction)) +
  geom_histogram(binwidth=2, position="stack", color="black", alpha=0.8) +
  geom_vline(xintercept=0, linetype="dashed", color="darkgray") +
  scale_fill_manual(values=c("Hypo in FB"="blue", "Hyper in FB"="red", "No change"="gray")) +
  theme_minimal(base_size=14) +
  labs(
    x = expression(Delta~Methylation~(Forebrain - Whole~brain)),
    y = "Number of regions",
    fill = "Direction",
    title = "Regional methylation differences"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face="bold")
  )



# Identify what regions are in tails, and go back to them in the regions file
# I added a "name" column that will allow us to identify it