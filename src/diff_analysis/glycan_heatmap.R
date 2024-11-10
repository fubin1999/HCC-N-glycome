library(tidyverse)
library(patchwork)

# Read data-----
# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# ancova_result <- read_csv("results/data/diff_analysis/glycan_ancova.csv")

abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
ancova_result <- read_csv(snakemake@input[[3]])

# Heatmap-----
diff_glycans <- ancova_result %>%
  filter(Effect == "group", p.adj < 0.05) %>%
  distinct(glycan) %>%
  pull(glycan)

heatmap_data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  filter(glycan %in% diff_glycans) %>%
  mutate(value = log(value + 1)) %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  group_by(glycan) %>%
  mutate(value = as.double(scale(value))) %>%
  group_by(glycan, group) %>%
  summarise(mean_z_score = mean(value), .groups = "drop")

heatmap <- ggplot(heatmap_data, aes(glycan, group)) +
  geom_tile(aes(fill = mean_z_score), color = "white", linewidth = 0.8) +
  scale_fill_gradient2(high = "#D26F32", low = "#275D87", mid = "white") +
  labs(x = "", y = "", fill = "z-score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1.1),
    legend.direction = "horizontal",
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# tgutil::ggpreview(width = 10, height = 2.6)
ggsave(snakemake@output[[1]], plot = heatmap, width = 10, height = 2.6)