library(tidyverse)
library(ggsignif)
library(patchwork)

# pred <- read_csv("results/data/cor_with_clinical/liver_function_model_pred.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

pred <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- pred %>%
  left_join(groups, by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

plot_boxplots <- function (data, .value, .target) {
  ggplot(data, aes(group, {{ .value }})) +
    geom_boxplot(aes(color = group)) +
    geom_jitter(aes(color = group), width = 0.25, alpha = 0.5) +
    scale_color_manual(values = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A")) +
    geom_signif(
      comparisons = list(c("HC", "CHB"), c("HC", "LC"), c("HC", "HCC")),
      map_signif_level = TRUE,
      step_increase = 0.1, vjust = 0
    ) +
    guides(color = FALSE) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
}

plot_true_boxplots <- function (data, .target) {
  data %>%
    plot_boxplots(true, .target) +
    labs(y = "True value", title = str_glue("True {.target}"))
}

plot_pred_boxplots <- function (data, .target) {
  data %>%
    plot_boxplots(pred, .target) +
    labs(y = "Predicted value", title = str_glue("Predicted {.target}"))
}

plot_df <- data %>%
  nest_by(target) %>%
  mutate(
    pred_plots = list(plot_pred_boxplots(data, target)),
    true_plots = list(plot_true_boxplots(data, target))
  ) %>%
  select(-data)

pred_p <- reduce(plot_df$pred_plots, `+`) + plot_layout(nrow = 1, axes = "collect_y")
true_p <- reduce(plot_df$true_plots, `+`) + plot_layout(nrow = 1, axes = "collect_y")
final_p <- pred_p / true_p
# tgutil::ggpreview(width = 10, height = 6.5)
ggsave(snakemake@output[[1]], final_p, width = 10, height = 6.5)
