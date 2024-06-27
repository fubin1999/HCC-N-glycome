library(tidyverse)

# Read data-----
# ancova_result <- read_csv("results/data/diff_analysis/glycan_ancova.csv")
# post_hoc_result <- read_csv("results/data/diff_analysis/glycan_post_hoc.csv")
# fold_change <- read_csv("results/data/diff_analysis/glycan_fold_change.csv")

ancova_result <- read_csv(snakemake@input[[1]])
post_hoc_result <- read_csv(snakemake@input[[2]])
fold_change <- read_csv(snakemake@input[[3]])

# Prepare data-----
data <- ancova_result %>%
  filter(Effect == "group", p.adj < 0.05) %>%
  select(glycan) %>%
  left_join(post_hoc_result, by = "glycan") %>%
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
  )

# Plot-----
p <- ggplot(data, aes(glycan, comparison)) +
  geom_point(aes(fill = logFC, size = logp, alpha = signif), color = "grey30", shape = 21) +
  guides(alpha = "none") +
  labs(x = "", y = "", color = expression(log[2]~FC), size = expression(-log[10]~p)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    legend.box = "horizontal"
  ) +
  coord_equal() +
  scale_alpha_manual(values = c(0, 1)) +
  scale_size_continuous(range = c(1, 5)) +
  scale_fill_gradient2(high = "#D26F32", low = "#275D87", mid = "white")
tgutil::ggpreview(width = 10, height = 4)

ggsave(snakemake@output[[1]], plot = p, width = 10, height = 4)