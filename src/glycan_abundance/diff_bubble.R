library(tidyverse)

# Read data-----
# post_hoc_result <- read_csv("results/data/glycan_abundance/posthoc_for_glycans.csv")
# fold_change <- read_csv("results/data/glycan_abundance/fold_change.csv")
# row_order <- read_csv("results/data/glycan_abundance/glycan_clusters.csv")

post_hoc_result <- read_csv(snakemake@input[["post_hoc"]])
fold_change <- read_csv(snakemake@input[["fold_change"]])
row_order <- read_csv(snakemake@input[["row_order"]])

# Prepare data-----
data <- post_hoc_result %>%
  select(glycan, group1, group2, p.adj) %>%
  inner_join(fold_change, by = c("glycan", "group1", "group2")) %>%
  mutate(
    comparison = str_glue("{group1} / {group2}"),
    comparison = factor(
      comparison,
      levels = c("HC / CHB", "HC / LC", "HC / HCC", "CHB / LC", "CHB / HCC", "LC / HCC")
    )
  ) %>%
  select(-group1, -group2) %>%
  mutate(
    logFC = log2(FC),
    logp = -log10(p.adj),
    signif = p.adj < 0.05,
    logp = if_else(is.infinite(logp), NA, logp),
    logp = if_else(is.na(logp), max(logp, na.rm = TRUE), logp)
  ) %>%
  left_join(row_order %>% mutate(rank = row_number()), by = "glycan")

# Plot-----
p <- ggplot(data, aes(comparison, reorder(glycan, desc(rank)))) +
  geom_point(aes(color = logFC, size = logp, alpha = signif)) +
  guides(alpha = "none") +
  labs(x = "", y = "", color = expression(log[2]~FC), size = expression(-log[10]~p)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    axis.text.y = element_blank(),
    legend.direction = "horizontal"
  ) +
  scale_alpha_manual(values = c(0, 1)) +
  scale_size_continuous(range = c(1, 3.5)) +
  scale_color_gradient2(high = "#CD0000", low = "#27408B", mid = "white")

# tgutil::ggpreview(width = 3.8, height = 7)
ggsave(snakemake@output[[1]], plot = p, width = 3.8, height = 7)