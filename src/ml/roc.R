source("renv/activate.R")

library(tidyverse)
library(pROC)
library(patchwork)

# Read data-----
# predictions <- read_csv("results/data/ml/predictions.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv") |>
#   select(sample, AFP)

predictions <- read_csv(snakemake@input[["predictions"]])
groups <- read_csv(snakemake@input[["groups"]])
clinical <- read_csv(snakemake@input[["clinical"]]) |>
  select(sample, AFP)

# Prepare data-----
temp_data <- predictions %>%
  select(-prediction) %>%
  left_join(groups, by = "sample") %>%
  left_join(clinical, by = "sample") %>%
  mutate(model = case_match(
    model, "HCC Fusion" ~ "fusion", "HCC Slim" ~ "slim"
  )) %>%
  pivot_wider(names_from = model, values_from = probability)

temp_data <- bind_rows(list(
  `Control/HCC` = temp_data %>% mutate(group = if_else(group == "HCC", "HCC", "Control")),
  `HC/HCC` = temp_data %>% filter(group %in% c("HC", "HCC")),
  `CHB/HCC` = temp_data %>% filter(group %in% c("CHB", "HCC")),
  `LC/HCC` = temp_data %>% filter(group %in% c("LC", "HCC"))
), .id = "comparison") %>%
  mutate(comparison = factor(comparison, levels = c("Control/HCC", "HC/HCC", "CHB/HCC", "LC/HCC")))

prepared_data <- bind_rows(list(
  `Test Set` = temp_data,
  `Test Set (AFP-)` = temp_data %>% filter(AFP < 10)
), .id = "subtype")

# Perform ROC analysis-----
roc_result <- prepared_data %>%
  group_by(subtype, comparison) %>%
  nest() %>%
  mutate(
    roc_fusion = map(data, ~ roc(data = .x, "target", "fusion", ci = TRUE)),
    roc_slim = map(data, ~ roc(data = .x, "target", "slim", ci = TRUE)),
    roc_AFP = map(data, ~ roc(data = .x, "target", "AFP", ci = TRUE))
  ) %>%
  ungroup() %>%
  select(-data) %>%
  pivot_longer(
     -c(subtype, comparison),
     names_to = "predictor",
     values_to = "roc",
     names_prefix = "roc_"
  ) %>%
  filter(!(subtype == "Test Set (AFP-)" & predictor == "AFP")) %>%
  mutate(predictor = case_match(
    predictor, "fusion" ~ "HCC Fusion", "slim" ~ "HCC Slim", "AFP" ~ "AFP"
  ))

# Tidy ROC AUC result-----
roc_auc <- roc_result  %>%
  mutate(
    auc = map_dbl(roc, ~ auc(.x)),
    ci_upper = map_dbl(roc, ~ .x$ci[3]),
    ci_lower = map_dbl(roc, ~ .x$ci[1])
  ) %>%
  select(-roc)

# Plot ROC curves-----
draw_roc_curve <- function (roc_list, .title) {
  ggroc(roc_list) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey") +
    labs(color = "", title = .title) +
    theme_bw() +
    theme(
      legend.position = c(0.7, 0.25),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.background = element_blank(),
      title = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#CC5F5A", "#7A848D", "#A2AFA6"))
}

roc_plots <- roc_result %>%
  nest_by(subtype, comparison) %>%
  mutate(
    roc_list = list(deframe(data)),
    plot = list(draw_roc_curve(roc_list, str_glue("{subtype}: {comparison}")))
  ) %>%
  select(-data, -roc_list) %>%
  ungroup()

final_plot <- reduce(roc_plots$plot, `+`) + plot_layout(nrow = 2)
# tgutil::ggpreview(width = 12, height = 6)

# Save result------
ggsave(snakemake@output[[1]], plot = final_plot, width = 12, height = 6)
write_csv(roc_auc, snakemake@output[[2]])