library(tidyverse)
library(patchwork)
library(ggprism)

# clinical <- read_csv("results/data/prepared/clinical.csv")
# subtypes <- read_csv("results/data/subtypes/consensus_cluster_result.csv")
# kw_result <- read_csv("results/data/subtypes/clinical_diff_kw.csv")
# post_hoc_result <- read_csv("results/data/subtypes/clinical_diff_post_hoc.csv")

clinical <- read_csv(snakemake@input[[1]])
subtypes <- read_csv(snakemake@input[[2]])
kw_result <- read_csv(snakemake@input[[3]])
post_hoc_result <- read_csv(snakemake@input[[4]])

plot_data <- clinical %>%
  select(-sex, -age) %>%
  pivot_longer(-sample, names_to = "clinical_variable", values_to = "value") %>%
  right_join(subtypes, by = "sample") %>%
  mutate(
    class = str_c("S", class),
    class = as.factor(class),
    value = log(value + 0.1)
  ) %>%
  semi_join(kw_result %>% filter(p.adj < 0.05), by = c("clinical_variable" = "feature"))

y_pos_df <- plot_data %>%
  group_by(clinical_variable) %>%
  summarise(y_max = max(value), y_range = max(value) - min(value))

p_value_data <- post_hoc_result %>%
  right_join(y_pos_df, by = c("feature" = "clinical_variable")) %>%
  group_by(feature) %>%
  mutate(y.position = y_max + (group2 - group1) * y_range * 0.08) %>%
  select(feature, group1, group2, label = p.adj.signif, y.position) %>%
  ungroup()

nested_plot_data <- plot_data %>%
  nest_by(clinical_variable)
nested_p_value_data <- p_value_data %>%
  nest_by(feature)

colors <- c("#1b9e77", "#d95f02", "#7570b3")

draw_boxplots <- function (value_data, p_data, .feature) {
  value_data %>%
    ggplot(aes(class, value)) +
    geom_boxplot(aes(color = class, fill = class), alpha = 0.1) +
    add_pvalue(p_data, bracket.shorten = 0.1, bracket.size = 0.4) +
    guides(color = "none", fill = "none") +
    labs(x = "Glycan Subtypes", y = str_glue("log({.feature})")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
}

plot_df <- nested_plot_data %>%
  left_join(
    nested_p_value_data,
    by = c("clinical_variable" = "feature"),
    suffix = c("_value", "_p")
  ) %>%
  mutate(plot = list(draw_boxplots(data_value, data_p, clinical_variable)))

p <- reduce(plot_df$plot, `+`) + plot_layout(nrow = 2)
# tgutil::ggpreview(p, width = 8, height = 5)
ggsave(snakemake@output[[1]], p, width = 8, height = 5)