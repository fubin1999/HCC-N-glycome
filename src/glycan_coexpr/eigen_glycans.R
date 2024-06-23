library(tidyverse)
library(broom)
library(rstatix)

abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clusters <- read_csv(snakemake@input[[3]]) %>%
  mutate(cluster = factor(as.integer(cluster)))

prepared_data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  right_join(clusters, by = "glycan") %>%
  group_by(glycan) %>%
  mutate(
    value = log2(value + 1),
    value = as.double(scale(value))
  ) %>%
  ungroup()

eigen_glycans <- prepared_data %>%
  nest_by(cluster) %>%
  mutate(
    data = list(pivot_wider(data, names_from = glycan, values_from = value)),
    data = list(column_to_rownames(data, "sample")),
    pca_fit = list(prcomp(data)),
    coords = list(augment(pca_fit))
  ) %>%
  select(cluster, coords) %>%
  unnest(coords) %>%
  select(sample = .rownames, cluster, eigen_glycan = .fittedPC1)

cluster_means <- prepared_data %>%
  summarise(mean = mean(value), .by = c(sample, cluster))

cor_eigen_with_mean <- cluster_means %>%
  left_join(eigen_glycans, by = c("sample", "cluster")) %>%
  group_by(cluster) %>%
  cor_test(mean, eigen_glycan, method = "spearman")

final_result <- eigen_glycans %>%
  left_join(cor_eigen_with_mean %>% select(cluster, cor)) %>%
  mutate(eigen_glycan = if_else(cor > 0, eigen_glycan, -eigen_glycan)) %>%
  select(-cor)

write_csv(final_result, snakemake@output[[1]])