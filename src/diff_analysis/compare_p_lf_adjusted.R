library(tidyverse)
library(cowplot)

# Load data-----
# no_adj_data <- read_csv("results/data/diff_analysis/glycan_ancova.csv")
# adj_data <- read_csv("results/data/diff_analysis/glycan_ancova_lf_adjusted.csv")

no_adj_data <- read_csv(snakemake@input[[1]])
adj_data <- read_csv(snakemake@input[[2]])

data <- bind_rows(list(adjusted = adj_data, not_adjusted = no_adj_data), .id = "adjustment") %>%
  filter(Effect == "group")

# Plot-----
scatter_p <- data %>%
  mutate(neglog10p = -log10(p.adj)) %>%
  select(adjustment, glycan, neglog10p) %>%
  pivot_wider(names_from = adjustment, values_from = neglog10p) %>%
  ggplot(aes(x = not_adjusted, y = adjusted)) +
  geom_point(shape = 21, color = "#275D87", size = 1.5) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    x = expression(-log[10]*p~"(raw)"),
    y = expression(-log[10]*p~"(adjusted)"),
  ) +
  theme_classic()

dot_plot <- data %>%
  mutate(neglog10p = -log10(p.adj)) %>%
  sjmisc::rec(adjustment, rec = "adjusted=ALBI-adjusted;not_adjusted=raw", suffix = "") %>%
  ggplot(aes(reorder(glycan, desc(neglog10p)), neglog10p)) +
  geom_line(aes(group = glycan), color = "grey90") +
  geom_point(aes(color = adjustment), size = 2, shape = 16) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("#D26F32", "#275D87")) +
  annotate("text", x = 1, y = -log10(0.05), label = "p = 0.05", hjust = 0, vjust = -1) +
  labs(x = "glycans", y = expression(-log[10]*p), color = NULL, fill = NULL) +
  guides(
    color = guide_legend(position = "inside"),
  ) +
  theme_classic() +
  theme(
    # axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.3, 0.9)
  )

# Put scatter plot onto dot plot-----
merged_p <- ggdraw(dot_plot) +
  draw_plot(scatter_p, x = 0.55, y = 0.35, width = 0.4, height = 0.6)
# tgutil::ggpreview(merged_p, width = 8, height = 4)
ggsave(snakemake@output[[1]], plot = merged_p, width = 8, height = 4)