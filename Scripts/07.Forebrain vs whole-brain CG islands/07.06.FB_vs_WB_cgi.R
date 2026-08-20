# This script compares methylation proportions of forebrain (ONT) and whole-brain (RRBS) datasets in CpG islands. 
# It computes Pearson and Spearman correlations, visualizes the relationship with scatter plots, 
# and analyzes what CpG islands are outliers based on prediction intervals from a linear regression model.

library(dplyr) # version 1.1.4
library(tidyr) # version 0.0.6


regions_name<- "cgi"

home <- path.expand("~")
path_fb <- paste0(home, "/Chaterjee/methylation_cgi/forebrain/")
path_wb <- paste0(home, "/Chaterjee/methylation_cgi/whole-brain/")

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
  # keep only rows where methylation and coverage exist
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
    df$meth <- as.numeric(df$meth) / 100 # To match forebrain dataset (percentage to fraction)
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

# 4. Compute Pearson and Spearman correlations
cor_pearson <- cor(merged_filtered$mean_meth_FB, merged_filtered$mean_meth_WB, method="pearson")
cor_spearman <- cor(merged_filtered$mean_meth_FB, merged_filtered$mean_meth_WB, method="spearman")

# 5. Fit linear regression and identify outliers
# Fit linear regression
lm_fit <- lm(mean_meth_FB ~ mean_meth_WB, data = merged_filtered)

# Add predicted values and 95% prediction intervals
pred <- predict(lm_fit,
                newdata = merged_filtered,
                interval = "prediction",
                level = 0.95)

merged_filtered <- merged_filtered %>%
  mutate(
    predicted = pred[, "fit"],
    PI_lower = pred[, "lwr"],
    PI_upper = pred[, "upr"],
    outside_PI = mean_meth_FB < PI_lower | mean_meth_FB > PI_upper
  )

table(merged_filtered$outside_PI)

# Identify outlying CpG islands based on the prediction intervals
outlier_cgis <- merged_filtered %>%
  filter(outside_PI) %>%
  arrange(desc(abs(mean_meth_FB - predicted)))

outlier_cgis
out_file <- paste0(home, "/Chaterjee/outlier_cgi_FB_vs_WB.txt")
write.table(outlier_cgis, out_file, sep = "\t", quote = F, col.names = T, row.names = F)

# Prediction grid
new_x <- data.frame(
  mean_meth_WB = seq(
    min(merged_filtered$mean_meth_WB),
    max(merged_filtered$mean_meth_WB),
    length.out = 200
  )
)

pred_grid <- predict(
  lm_fit,
  newdata = new_x,
  interval = "prediction",
  level = 0.95
)

pred_grid <- cbind(new_x, pred_grid)

pred_grid <- pred_grid %>%
  mutate(
    lwr = pmax(lwr, 0),
    upr = pmin(upr, 1)
  )

# Plot
p <- ggplot(merged_filtered,
            aes(x = mean_meth_WB, y = mean_meth_FB)) +
  
  geom_ribbon(
    data = pred_grid,
    aes(x = mean_meth_WB, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.10
  ) +
  
  geom_line(
    data = pred_grid,
    aes(x = mean_meth_WB, y = fit),
    inherit.aes = FALSE,
    color = "black"
  ) +
  
  geom_point(
    aes(color = outside_PI, alpha = outside_PI)
  ) +
  
  scale_color_manual(
    values = c("FALSE" = "black",
               "TRUE" = "#1f4e99"),
    labels = c("FALSE" = "Within prediction interval",
               "TRUE" = "Outside prediction interval"),
    name = NULL
  ) +
  
  scale_alpha_manual(
    values = c("FALSE" = 0.1,
               "TRUE" = 0.9),
    guide = "none"
  ) +
  
  xlab("Whole-brain RRBS") +
  ylab("Forebrain ONT") +
  theme_classic() +
  
  theme(legend.position = "none")

p
out_file <- paste0(home, "/Chaterjee/cgi_fb_vs_wb.pdf")
ggsave(out_file, p)