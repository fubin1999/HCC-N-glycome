source("renv/activate.R")

library(tidyverse)
library(patchwork)

# trait_data <- read_csv("results/data/derived_traits/filtered_derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# post_hoc_result <- read_csv("results/data/derived_traits/posthoc_for_derived_traits.csv")
# fold_change <- read_csv("results/data/derived_traits/fold_change.csv")

trait_data <- read_csv(snakemake@input[["traits"]])
groups <- read_csv(snakemake@input[["groups"]])
post_hoc_result <- read_csv(snakemake@input[["post_hoc"]])
fold_change <- read_csv(snakemake@input[["fold_change"]])

bubble_data <- fold_change %>%
  inner_join(post_hoc_result, by = c("trait", "group1", "group2")) %>%
  mutate(comparison = paste0(group1, " / ", group2)) %>%
  mutate(signif = if_else(p.adj < 0.05, TRUE, FALSE)) %>%
  mutate(logFC = log2(FC), logp = -log10(p.adj)) %>%
  mutate(
    logp = if_else(is.infinite(logp), NA, logp),
    logp = if_else(is.na(logp), max(logp, na.rm = TRUE), logp)
  ) %>%
  select(trait, comparison, logFC, logp, signif)

bubble_p <- ggplot(bubble_data, aes(trait, comparison)) +
  geom_point(aes(fill = logFC, size = logp, alpha = signif), color = "grey30", shape = 21) +
  guides(alpha = "none") +
  labs(x = "", y = "", color = expression(log[2]~FC), size = expression(-log[10]~p)) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    legend.box = "horizontal"
  ) +
  scale_alpha_manual(values = c(0, 1)) +
  scale_size_continuous(range = c(1, 3.5)) +
  scale_fill_gradient2(high = "#CD0000", low = "#27408B", mid = "white")

heatmap_data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  semi_join(post_hoc_result, by = "trait") %>%
  group_by(trait) %>%
  mutate(value = as.vector(scale(value))) %>%
  ungroup()

heatmap_p <- ggplot(heatmap_data, aes(trait, group, fill = value)) +
  geom_tile() +
  labs(x = "", y = "", fill = "Z-score") +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gradient2(
    high = "#CD0000", low = "#27408B", mid = "white",
    breaks = c(-2, 0, 2), limits = c(-2, 2)
  )

heatmap_p / bubble_p + plot_layout(guides = "collect") &
  theme(legend.direction = "horizontal", legend.position = "right")
# tgutil::ggpreview(width = 12, height = 4)
ggsave(snakemake@output[[1]], width = 12, height = 4)