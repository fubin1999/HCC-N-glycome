library(tidyverse)
library(yardstick)
library(patchwork)
library(gt)


# Read data-----
preds <- read_csv("results/data/ml/preds.csv")
preds <- preds %>% 
  mutate(true = factor(true), pred = factor(pred)) %>% 
  mutate(model = case_match(
    model, 
    "HCC_HC" ~ "HCC vs HC",
    "HCC_CHB" ~ "HCC vs CHB",
    "HCC_LC" ~ "HCC vs LC",
    "global" ~ "HCC vs Rest",
  )) %>% 
  mutate(model = factor(model, levels = c("HCC vs HC", "HCC vs CHB", "HCC vs LC", "HCC vs Rest")))

# Metrics-----
class_metrics <- metric_set(accuracy, bal_accuracy, f_meas, sensitivity, specificity)
metrics1 <- preds %>% 
  group_by(model, dataset) %>% 
  class_metrics(true, estimate = pred)

metrics2 <- preds %>% 
  group_by(model, dataset) %>% 
  roc_auc(true, proba, event_level = "second")

metrics <- bind_rows(metrics1, metrics2) %>% 
  select(model, dataset, metric = .metric, score = .estimate)

make_table <- function(data) {
  data %>% 
    pivot_wider(names_from = metric, values_from = score) %>% 
    gt(rowname_col = "dataset") %>% 
    fmt_number() %>% 
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_title()
    ) %>%
    opt_stylize(style = 1)
}

table_df <- metrics %>% 
  mutate(
    metric = case_match(
      metric,
      "accuracy" ~ "ACC",
      "bal_accuracy" ~ "BACC",
      "f_meas" ~ "F1",
      "sensitivity" ~ "SEN",
      "specificity" ~ "SPE",
      "roc_auc" ~ "AUC",
    ),
    dataset = case_match(
      dataset,
      "test" ~ "Int. Val.",
      "val1" ~ "Ext. Val. 1",
      "val2" ~ "Ext. Val. 2"
    )
  ) %>% 
  nest_by(model) %>% 
  mutate(table = list(make_table(data) %>% tab_header(title = model))) %>% 
  select(-data) %>% 
  mutate(filename = paste0(
    "results/figures/ml/metrics_", 
    janitor::make_clean_names(model), 
    ".html")
  )

walk2(table_df$table, table_df$filename, ~gtsave(.x, filename = .y))

roc_auc_ci <- preds %>% 
  nest_by(model, dataset) %>% 
  mutate(
    roc = list(pROC::roc(data$true, data$proba, ci = TRUE, ci.method = "delong")),
    auc = roc$ci[[2]],
    ci_lower = roc$ci[[1]],
    ci_upper = roc$ci[[3]]
  ) %>% 
  select(model, dataset, auc, ci_lower, ci_upper)

write_csv(roc_auc_ci, "results/data/ml/roc_auc_ci.csv")

# Plot ROC Curves-----
plot_roc <- function(data) {
  data %>% 
    group_by(dataset) %>% 
    roc_curve(true, proba, event_level = "second") %>% 
    autoplot() +
    guides(color = guide_legend(title = NULL, position = "inside")) +
    scale_color_manual(
      values = c(val1 = "#3A3D8F", val2 = "#5DA5DA", test = "#ED6A5A"),
      labels = c(test = "Int. Val.", val1 = "Ext. Val. 1", val2 = "Ext. Val. 2"),
    ) +
    theme(
      panel.grid = element_blank(),
      legend.position.inside = c(0.7, 0.25),
      plot.title = element_text(hjust = 0.5)
    )
}

roc_df <- preds %>% 
  nest_by(model) %>% 
  mutate(plot = list(plot_roc(data) + ggtitle(model))) %>% 
  arrange(model)
combined_roc <- reduce(roc_df$plot, `+`) + plot_layout(nrow = 1)

tgutil::ggpreview(combined_roc, width = 12, height = 3)
ggsave("results/figures/ml/roc_curves.pdf", combined_roc, width = 12, height = 3)

# Plot PR Curves-----
plot_pr_curve <- function(data) {
  data %>% 
    group_by(dataset) %>%
    pr_curve(true, proba, event_level = "second") %>% 
    autoplot() +
    guides(color = guide_legend(title = NULL, position = "inside")) +
    scale_color_manual(
      values = c(val1 = "#3A3D8F", val2 = "#5DA5DA", test = "#ED6A5A"),
      labels = c(test = "Int. Val.", val1 = "Ext. Val. 1", val2 = "Ext. Val. 2"),
    ) +
    theme(
      panel.grid = element_blank(),
      legend.position.inside = c(0.3, 0.25),
      plot.title = element_text(hjust = 0.5)
    )
}

pr_df <- preds %>% 
  nest_by(model) %>% 
  mutate(plot = list(plot_pr_curve(data) + ggtitle(model))) %>% 
  arrange(model)
combined_pr <- reduce(pr_df$plot, `+`) + plot_layout(nrow = 1)

tgutil::ggpreview(combined_pr, width = 12, height = 3)
ggsave("results/figures/ml/pr_curves.pdf", combined_pr, width = 12, height = 3)
