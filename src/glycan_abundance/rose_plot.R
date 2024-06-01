library(tidyverse)

# post_hoc_result <- read_csv("results/data/glycan_abundance/posthoc_for_glycans.csv")
post_hoc_result <- read_csv(snakemake@input[[1]])

plot_data <- post_hoc_result %>%
  mutate(
    comparison = str_glue("{group1} / {group2}"),
    signif = p.adj < 0.05,
    regulate = if_else(estimate > 0, "up", "down")
  ) %>%
  select(glycan, comparison, signif, regulate) %>%
  summarise(n = sum(signif), .by = c(comparison, regulate))

ggplot(plot_data, aes(comparison, n)) +
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
ggsave(snakemake@output[[1]], width = 3, height = 3.5)