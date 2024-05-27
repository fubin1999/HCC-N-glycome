source("renv/activate.R")

library(tidyverse)

# post_hoc_result <- read_csv("results/data/derived_traits/posthoc_for_derived_traits.csv")
# fold_change <- read_csv("results/data/derived_traits/fold_change.csv")
post_hoc_result <- read_csv(snakemake@input[[1]])
fold_change <- read_csv(snakemake@input[[2]])

data <- fold_change %>%
  inner_join(post_hoc_result, by = c("trait", "group1", "group2")) %>%
  mutate(comparison = paste0(group1, " / ", group2)) %>%
  mutate(signif = if_else(p.adj < 0.05, TRUE, FALSE)) %>%
  mutate(logFC = log2(FC), logp = -log10(p.adj)) %>%
  select(trait, comparison, logFC, logp, signif)

ggplot(data, aes(trait, comparison)) +
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
tgutil::ggpreview(width = 10, height = 2.5)
ggsave(snakemake@output[[1]], width = 10, height = 2.5)