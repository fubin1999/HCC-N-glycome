library(tidyverse)
library(broom)

abundance <- read_csv("results/data/prepared/processed_abundance.csv")
clusters <- read_csv("results/data/subtypes/consensus_cluster_result.csv")

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
  select(sample, subtype = cluster, PC1 = .fittedPC1, PC2 = .fittedPC2) %>% 
  mutate(subtype = paste0("S", subtype))

p <- ggplot(plot_data, aes(PC1, PC2)) +
  geom_point(aes(fill = subtype), size = 2, shape = 21, color = "black") +
  labs(fill = NULL) +
  scale_fill_manual(values = c(S1 = "#1b9e77", S2 = "#d95f02", S3 = "#7570b3")) +
  guides(fill = guide_legend(position = "inside")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position.inside = c(0.16, 0.78)
  )
tgutil::ggpreview(plot = p, width = 3, height = 3)

ggsave("results/figures/subtypes/subtype_pca.pdf", plot = p, width = 3, height = 3)
