library(tidyverse)

ancova_result <- read_csv("results/data/diff_analysis/glycan_ancova.csv")
post_hoc_result <- read_csv("results/data/diff_analysis/glycan_post_hoc.csv")

ancova_result <- read_csv(snakemake@input[[1]])
post_hoc_result <- read_csv(snakemake@input[[2]])

diff_glycans <- ancova_result %>%
  filter(Effect == "group", p.adj < 0.05) %>%
  pull(glycan)

plot_data <- post_hoc_result %>%
  filter(glycan %in% diff_glycans) %>%
  mutate(
    comparison = str_glue("{group1} / {group2}"),
    signif = p.adj < 0.05,
    regulate = if_else(estimate > 0, "up", "down")
  ) %>%
  select(glycan, comparison, signif, regulate) %>%
  summarise(n = sum(signif), .by = c(comparison, regulate))

p <- ggplot(plot_data, aes(comparison, n)) +
  geom_col(aes(fill = regulate), alpha = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), color = "white") +
  coord_polar() +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    #panel.grid.major.y = element_blank(),
    axis.text.x = element_text(vjust = -1),
    axis.text.y = element_blank(),
    legend.position = "bottom",
  ) +
  scale_fill_manual(values = c("up" = "#CC5F5A", "down" = "#97A5C0"))
# tgutil::ggpreview(width = 4, height = 5)
ggsave(snakemake@output[[1]], plot = p, width = 3, height = 3.5)