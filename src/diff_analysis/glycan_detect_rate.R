library(tidyverse)
library(rstatix)
library(ggprism)
library(patchwork)

# raw_abundance <- read_csv("results/data/prepared/raw_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
raw_abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

detection <- raw_abundance %>%
  right_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  pivot_longer(-c(sample, group), names_to = "glycan", values_to = "value") %>%
  mutate(detected = !is.na(value)) %>%
  group_by(group, glycan, detected) %>%
  count() %>%
  ungroup() %>%
  complete(group, glycan, detected, fill = list(n = 0))

tidy_fisher_test <- function (data) {
  data %>%
    pivot_wider(names_from = group, values_from = n) %>%
    column_to_rownames("detected") %>%
    as.matrix() %>% as.table() %>%
    pairwise_fisher_test(detailed = TRUE, p.adjust.method = "BH") %>%
    as_tibble() %>%
    select(group1, group2, estimate, conf.low, conf.high, p)
}

fisher_result <- detection %>%
  nest_by(glycan) %>%
  mutate(fisher_result = list(tidy_fisher_test(data))) %>%
  select(-data) %>%
  unnest(fisher_result) %>%
  adjust_pvalue()

write_csv(fisher_result, snakemake@output[[1]])

diff_glycans <- fisher_result %>%
  filter(p.adj < 0.05) %>%
  distinct(glycan) %>%
  pull(glycan)

plot_data <- detection %>%
  filter(glycan %in% diff_glycans) %>%
  group_by(glycan, group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    detected = if_else(detected, "Yes", "No"),
    detected = factor(detected, levels = c("Yes", "No")),
    n_label_mid = if_else(prop > 0.08, n, NA),
    n_label_top = if_else(is.na(n_label_mid), n, NA)
  )

step_increase <- 0.8
p_val_data <- fisher_result %>%
  add_significance("p.adj") %>%
  mutate(
    label = p.adj.signif,
    y.position = case_when(
      (group1 == "HC") & (group2 == "CHB") ~ 1.0 + step_increase * 0.1,
      (group1 == "CHB") & (group2 == "LC") ~ 1.0 + step_increase * 0.1,
      (group1 == "LC") & (group2 == "HCC") ~ 1.0 + step_increase * 0.1,
      (group1 == "HC") & (group2 == "LC") ~ 1.0 + step_increase * 0.2,
      (group1 == "CHB") & (group2 == "HCC") ~ 1.0 + step_increase * 0.35,
      (group1 == "HC") & (group2 == "HCC") ~ 1.0 + step_increase * 0.45
    )
  ) %>%
  select(glycan, group1, group2, label, y.position)

plot_fisher <- function (data, p_val_data, .title) {
  data %>%
    ggplot(aes(group, prop)) +
    geom_col(aes(fill = detected)) +
    scale_fill_manual(values = c("Yes" = "#E64B35FF", "No" = "#4DBBD5FF")) +
    geom_text(
      aes(label = n_label_mid),
      position = position_stack(vjust = 0.5),
      color = "white"
    ) +
    geom_text(
      aes(label = n_label_top),
      #position = position_stack(vjust = 1),
      nudge_y = 0.05,
      color = "white"
    ) +
    add_pvalue(data = p_val_data, bracket.size = 0.3, bracket.shorten = 0.05) +
    labs(
      fill = "Detected",
      title = .title
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

plot_df <- plot_data %>%
  nest_by(glycan, .key = "data") %>%
  left_join(p_val_data %>% nest_by(glycan, .key = "p_val_data")) %>%
  mutate(plot = list(plot_fisher(data, p_val_data, glycan))) %>%
  select(glycan, plot)

final_p <- reduce(plot_df$plot, `+`) +
  plot_layout(nrow = 2, guides = "collect") &
  theme(
    legend.position = "right"
  )
tgutil::ggpreview(plot = final_p, width = 14, height = 7)
ggsave(snakemake@output[[2]], plot = final_p, width = 14, height = 7)