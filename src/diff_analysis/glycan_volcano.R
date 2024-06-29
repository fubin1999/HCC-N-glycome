library(tidyverse)
library(ggrepel)
library(patchwork)


ancova_result <- read_csv("results/data/diff_analysis/glycan_ancova.csv")
post_hoc_result <- read_csv("results/data/diff_analysis/glycan_post_hoc.csv")
fold_change <- read_csv("results/data/diff_analysis/glycan_fold_change.csv")

ancova_result <- read_csv(snakemake@input[[1]])
post_hoc_result <- read_csv(snakemake@input[[2]])
fold_change <- read_csv(snakemake@input[[3]])

diff_glycans <- ancova_result %>%
  filter(Effect == "group", p.adj < 0.05) %>%
  pull(glycan)

data <- post_hoc_result %>%
  mutate(ancova_signif = glycan %in% diff_glycans) %>%
  left_join(fold_change, by = c("glycan", "group1", "group2")) %>%
  select(glycan, group1, group2, ancova_signif, p.adj, FC)

plot_data <- data %>%
  mutate(logp = -log10(p.adj), logFC = log2(FC)) %>%
  mutate(regulate = case_when(
    ancova_signif & p.adj < 0.05 & logFC > log2(1) ~ "up",
    ancova_signif & p.adj < 0.05 & logFC < -log2(1) ~ "down",
    .default = "no"
  )) %>%
  mutate(regulate = factor(regulate, levels = c("up", "down", "no"))) %>%
  filter(group1 == "HC") %>%
  mutate(
    comparison = str_c(group2, " vs ", group1),
    comparison = factor(comparison, levels = c("CHB vs HC", "LC vs HC", "HCC vs HC"))
  ) %>%
  select(glycan, comparison, logp, logFC, regulate)

ylim <- c(0, max(plot_data$logp))

plot_volcano <- function (data, .title) {
  ggplot(data, aes(logFC, logp)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_point(aes(color = regulate, size = abs(logFC)), alpha = 0.5) +
    ggtitle(.title)
}

plot_df <- plot_data %>%
  nest_by(comparison) %>%
  mutate(plot = list(plot_volcano(data, comparison))) %>%
  select(-data)

p <- reduce(plot_df$plot, `+`) +
  plot_layout(nrow = 1, guides = "collect") &
  lims(x = c(-0.8, 1.4), y = ylim) &
  scale_color_manual(
    values = c("up" = "#CD0000", "down" = "steelblue", "no" = "grey"),
  ) &
  scale_size_continuous(range = c(1, 5), limits = c(0, 1.5)) &
  labs(x = "log2FC", y = "-log10p", color = "Regulate", size = "|log2FC|") &
  theme_classic()
# tgutil::ggpreview(plot = p, width = 10, height = 3)

ggsave(snakemake@output[[1]], plot = p, width = 10, height = 3)