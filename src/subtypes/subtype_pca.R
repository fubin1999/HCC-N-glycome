library(tidyverse)
library(broom)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# clusters <- read_csv("results/data/subtypes/consensus_cluster_result.csv")

abundance <- read_csv(snakemake@input[[1]])
clusters <- read_csv(snakemake@input[[2]])

data <- clusters %>%
  mutate(cluster = factor(class), .keep = "unused") %>%
  left_join(abundance, by = "sample")

pca_fit <- data %>%
  select(-sample, -cluster) %>%
  prcomp(scale = TRUE)

plot_data <- pca_fit %>%
  augment(data) %>%
  select(sample, cluster, PC1 = .fittedPC1, PC2 = .fittedPC2)

p <- ggplot(plot_data, aes(PC1, PC2)) +
  geom_point(aes(color = cluster), alpha = 0.5) +
  stat_ellipse(aes(fill = cluster), geom = "polygon", alpha = 0.1) +
  scale_fill_manual(values = c(`1` = "#1b9e77", `2` = "#d95f02", `3` = "#7570b3")) +
  scale_color_manual(values = c(`1` = "#1b9e77", `2` = "#d95f02", `3` = "#7570b3")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
  )
tgutil::ggpreview(plot = p, width = 5, height = 4)

ggsave(snakemake@output[[1]], plot = p, width = 5, height = 4)