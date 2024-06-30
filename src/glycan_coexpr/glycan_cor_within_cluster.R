library(tidyverse)
library(rstatix)
library(ggsignif)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# clusters <- read_csv("results/data/glycan_coexpr/glycan_clusters.csv")

abundance <- read_csv(snakemake@input[[1]])
clusters <- read_csv(snakemake@input[[2]])

data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  mutate(value = log2(value)) %>%
  right_join(clusters, by = "glycan") %>%
  mutate(cluster = factor(cluster))

cor_result <- data %>%
  select(-cluster) %>%
  pivot_wider(names_from = glycan, values_from = value) %>%
  column_to_rownames("sample") %>%
  cor_mat() %>%
  rename(glycan1 = rowname) %>%
  pivot_longer(-glycan1, names_to = "glycan2", values_to = "r") %>%
  filter(glycan1 != glycan2)

cor_result_with_type <- cor_result %>%
  inner_join(clusters, by = c("glycan1" = "glycan")) %>%
  rename(cluster1 = cluster) %>%
  inner_join(clusters, by = c("glycan2" = "glycan")) %>%
  rename(cluster2 = cluster) %>%
  mutate(type = if_else(cluster1 == cluster2, "intra", "inter"))

p <- ggplot(cor_result_with_type, aes(type, r)) +
  geom_boxplot(color = "steelblue", fill = "steelblue", alpha = 0.1) +
  geom_signif(
    comparisons = list(c("inter", "intra")),
    map_signif_level = function(p) sprintf("p = %.2g", p),
    vjust = -0.2, tip_length = 0
  ) +
  labs(x = "Type of Correlation", y = "Pearson's r") +
  scale_x_discrete(labels = c("inter" = "Inter-GCM", "intra" = "Intra-GCM")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )
# tgutil::ggpreview(plot = p, width = 2.5, height = 3.5)

ggsave(snakemake@output[[1]], width = 2.5, height = 3.5)