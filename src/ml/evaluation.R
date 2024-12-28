library(tidyverse)
library(yardstick)
library(patchwork)
library(gt)


# Read data-----
preds <- read_csv("results/data/ml/preds.csv")
clinical <- read_csv("results/data/prepared/unfiltered_clinical.csv")

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


# Plot ROC comparing models and AFP-----
plot_roc_compared_with_AFP <- function(data, AFP_auc, model_auc) {
  data %>% 
    select(sample, true, model = proba, AFP) %>% 
    pivot_longer(c(model, AFP), names_to = "predictor", values_to = "value") %>% 
    group_by(predictor) %>% 
    roc_curve(true, value, event_level = "second") %>% 
    autoplot() +
    scale_color_manual(
      values = c(AFP = "#3A3D8F", model = "#ED6A5A"),
      labels = c(
        AFP = str_glue("AFP ({sprintf('%.3f', AFP_auc)})"), 
        model = str_glue("Model ({sprintf('%.3f', model_auc)})")
      )
    ) +
    guides(color = guide_legend(title = NULL, position = "inside")) +
    theme(
      panel.grid = element_blank(),
      legend.position.inside = c(0.65, 0.2),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

prepared_for_delong <- preds %>% 
  filter(dataset == "test") %>%
  select(-dataset) %>% 
  left_join(clinical %>% select(sample, AFP), by = "sample")

delong_res <- prepared_for_delong %>% 
  nest_by(model) %>% 
  mutate(
    model_roc = list(pROC::roc(data$true, data$proba)),
    AFP_roc = list(pROC::roc(data$true, data$AFP)),
    roc_test_res = list(pROC::roc.test(model_roc, AFP_roc)),
    delong_p = roc_test_res$p.value
  ) %>% 
  select(model, delong_p)

auc_df <- prepared_for_delong %>% 
  select(sample, true, proba, AFP, model) %>% 
  pivot_longer(c(proba, AFP), names_to = "predictor", values_to = "value") %>%
  group_by(model, predictor) %>% 
  roc_auc(true, value, event_level = "second") %>% 
  pivot_wider(names_from = predictor, values_from = .estimate, names_prefix = "auc_") %>% 
  select(-c(.metric, .estimator))

AFP_roc_plot_df <- prepared_for_delong %>% 
  nest_by(model) %>% 
  left_join(delong_res, by = "model") %>%
  left_join(auc_df, by = "model") %>%
  arrange(model) %>%
  mutate(plot = list(
    plot_roc_compared_with_AFP(data, auc_AFP, auc_proba) + 
      ggtitle(model) +
      labs(subtitle = paste("DeLong p-value:", round(delong_p, 3)))
  ))


AFP_ROC_p <- reduce(AFP_roc_plot_df$plot, `+`) + plot_layout(nrow = 1)
tgutil::ggpreview(AFP_ROC_p, width = 12, height = 4)
ggsave("results/figures/ml/roc_AFP.pdf", AFP_ROC_p, width = 12, height = 4)
