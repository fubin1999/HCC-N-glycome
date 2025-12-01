library(tidyverse)

# Load data-----
no_adj_data <- read_csv("results/data/diff_analysis/glycan_ancova.csv")
adj_data <- read_csv("results/data/diff_analysis/glycan_ancova_lf_adjusted.csv")

data <- bind_rows(list(adjusted = adj_data, not_adjusted = no_adj_data), .id = "adjustment") %>%
  filter(Effect == "group")

# Plot-----
glycans_to_plot <- data %>%
  select(adjustment, glycan, p.adj) %>%
  pivot_wider(names_from = adjustment, values_from = p.adj) %>%
  filter(adjusted < 0.05, not_adjusted < 0.05) %>%
  pull(glycan)

plot_data <- data %>%
  filter(Effect == "group") %>%
  mutate(neglog10p = -log10(p.adj)) %>%
  mutate(adjustment = factor(adjustment, levels = c("adjusted", "not_adjusted"))) %>% 
  mutate(label = if_else(neglog10p > -log10(0.05), glycan, ""))

dot_plot <- plot_data %>% 
  ggplot(aes(reorder(glycan, desc(neglog10p)), neglog10p)) +
  geom_line(aes(group = glycan), color = "grey90", size = 1.5) +
  geom_point(aes(color = adjustment), size = 3.5, shape = 16) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(
    values = c("#D26F32", "#275D87"),
    labels = c("adjusted" = "Adjusted for ALBI", "not_adjusted" = "Not adjusted for ALBI")
  ) +
  annotate("text", x = 1, y = -log10(0.05), label = "p = 0.05", hjust = 0, vjust = -0.6) +
  labs(x = "glycans", y = expression(-log[10]*p), color = NULL, fill = NULL) +
  guides(
    color = guide_legend(position = "inside"),
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.7, 0.8),
    legend.direction = "horizontal",
  )

# Put scatter plot onto dot plot-----
# tgutil::ggpreview(dot_plot, width = 6, height = 3)
ggsave(snakemake@output[[1]], plot = dot_plot, width = 6, height = 3)

write_csv(plot_data, "results/source_data/Supplementary_Figure_10.csv")
