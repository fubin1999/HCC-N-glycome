library(tidyverse)
library(patchwork)

# pred <- read_csv("results/data/cor_with_clinical/liver_function_model_pred.csv")
pred <- read_csv(snakemake@input[[1]])

pred_plot_data <- pred %>%
  select(-true) %>%
  group_by(target) %>%
  mutate(pred = as.double(scale(pred))) %>%
  ungroup() %>%
  pivot_wider(names_from = target, values_from = pred) %>%
  mutate(`-ALB` = -ALB, .keep = "unused") %>%
  rowwise() %>%
  mutate(color = sum(c_across(-sample))) %>%
  ungroup() %>%
  pivot_longer(-c(sample, color), names_to = "target", values_to = "pred")

true_plot_data <- pred %>%
  select(-pred) %>%
  group_by(target) %>%
  mutate(true = as.double(scale(true))) %>%
  ungroup() %>%
  pivot_wider(names_from = target, values_from = true) %>%
  mutate(`-ALB` = -ALB, .keep = "unused") %>%
  rowwise() %>%
  mutate(color = sum(c_across(-sample))) %>%
  ungroup() %>%
  pivot_longer(-c(sample, color), names_to = "target", values_to = "true")

# Remove outlier samples
outlier_samples <- true_plot_data %>%
  pivot_longer(-c(sample, target, color), names_to = "value_type", values_to = "value") %>%
  mutate(outlier = abs(value) > 5) %>%
  select(sample, outlier, target) %>%
  pivot_wider(names_from = target, values_from = outlier) %>%
  rowwise() %>%
  mutate(outlier = sum(c_across(-sample)) > 0) %>%
  ungroup() %>%
  filter(outlier) %>%
  pull(sample)
true_plot_data <- true_plot_data %>%
  filter(!sample %in% outlier_samples)

plot_para_coord <- function (data, .value_type) {
  ggplot(data, aes(target, {{ .value_type }}, group = sample)) +
    geom_line(aes(color = color), alpha = 0.5) +
    scale_color_gradient2(high = "#D26F32", mid = "white", low = "#275D87") +
    labs(x = "Clinical Variables", y = "Z-score") +
    guides(color = FALSE) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_line(color = "grey"),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)))
}

pred_p <- plot_para_coord(pred_plot_data, pred) +
  ggtitle("Predicted values")
true_p <- plot_para_coord(true_plot_data, true) +
  ggtitle("True values")

p <- pred_p + true_p
# tgutil::ggpreview(plot = p, width = 6, height = 3)
ggsave(snakemake@output[[1]], plot = p, width = 6, height = 3)