library(tidyverse)
library(patchwork)

post_hoc_result <- read_csv(snakemake@input[[1]])

bubble_data <- post_hoc_result %>%
  mutate(
    comparison = paste0(group1, " / ", group2),
    comparison = factor(comparison, levels = c(
      "HC / CHB", "HC / LC", "HC / HCC",
      "CHB / LC", "CHB / HCC", "LC / HCC"
    ))
  ) %>%
  mutate(signif = if_else(p.adj < 0.05, TRUE, FALSE)) %>%
  mutate(
    logp = -log10(p.adj),
    logp = if_else(is.infinite(logp), NA, logp),
    logp = if_else(is.na(logp), max(logp, na.rm = TRUE), logp)
  ) %>%
  select(trait, comparison, cohens_d, logp, signif)

bubble_p <- ggplot(bubble_data, aes(trait, comparison)) +
  geom_point(aes(fill = cohens_d, size = logp, alpha = signif), color = "grey30", shape = 21) +
  guides(alpha = "none") +
  labs(x = "", y = "", fill = "Cohen's d", size = expression(-log[10]~p)) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    axis.text.y = element_text(hjust = 1),
    legend.box = "horizontal"
  ) +
  scale_alpha_manual(values = c(0, 1)) +
  scale_size_continuous(range = c(1, 5)) +
  scale_fill_gradient2(high = "#D26F32", low = "#275D87", mid = "white") +
  scale_y_discrete(position = "right")

# tgutil::ggpreview(width = 10, height = 4)
ggsave(snakemake@output[[1]], plot = bubble_p, width = 10, height = 4)