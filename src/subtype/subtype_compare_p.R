library(tidyverse)
library(cowplot)

# ancova_result <- read_csv("results/data/subtype/subtype_ancova_with_HC.csv")
ancova_result <- read_csv(snakemake@input[[1]])

plot_data <- ancova_result %>%
  mutate(log_p = -log10(p.adj)) %>%
  select(glycan, group2, log_p) %>%
  pivot_wider(names_from = group2, values_from = log_p) %>%
  mutate(more_signif = if_else(HCC_S1 > HCC_S2, "HCC_S1", "HCC_S2"))

scatter_p <- ggplot(plot_data, aes(HCC_S2, HCC_S1)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(aes(color = more_signif), shape = 16) +
  scale_color_manual(values = c("HCC_S1" = "#E64B35FF", "HCC_S2" = "#4DBBD5FF")) +
  coord_equal() +
  lims(x = c(0, 17), y = c(0, 17)) +
  labs(
    x = "-log10(p), HCC_S2 vs HC",
    y = "-log10(p), HCC_S1 vs HC",
    color = "More Signif."
  ) +
  theme_classic()
# tgutil::ggpreview(plot = scatter_p, width = 4, height = 3)

bar_p <- plot_data %>%
  mutate(more_signif = str_sub(more_signif, 5, -1)) %>%
  group_by(more_signif) %>%
  count() %>%
  ggplot(aes(more_signif, n)) +
  geom_col(aes(fill = more_signif)) +
  geom_text(aes(label = n), nudge_y = 4.5) +
  scale_fill_manual(values = c("S1" = "#E64B35FF", "S2" = "#4DBBD5FF")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  guides(fill = "none") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.line.x = element_line()
  )
# tgutil::ggpreview(plot = bar_p, width = 1.5, height = 2.5)

final_p <- ggdraw(scatter_p) +
  draw_plot(bar_p, 0.4, 0.15, 0.3, 0.4)
# tgutil::ggpreview(plot = final_p, width = 4, height = 3)
ggsave(snakemake@output[[1]], plot = final_p, width = 4, height = 3)