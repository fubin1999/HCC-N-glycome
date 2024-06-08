library(tidyverse)
library(yardstick)
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
), .id = "subtype") %>%
  mutate(target = factor(target, levels = c(TRUE, FALSE))) %>%
  pivot_longer(c(AFP, fusion, slim), names_to = "predictor", values_to = "value") %>%
  mutate(predictor = case_match(
    predictor,
    "fusion" ~ "HCC Fusion",
    "slim" ~ "HCC Slim",
    "AFP" ~ "AFP"
  )) %>%
  mutate(predictor = factor(predictor, levels = c("HCC Fusion", "HCC Slim", "AFP"))) %>%
  filter(!((subtype == "Test Set (AFP-)") & (predictor == "AFP")))

# ROC Curves-----
plot_roc_curve <- function (data, .title) {
  data %>%
    roc_curve(target, value) %>%
    ggplot(aes(1 - specificity, sensitivity, color = predictor)) +
    geom_path() +
    geom_abline(lty = 3) +
    coord_equal() +
    ggtitle(.title) +
    guides(color = guide_legend(position = "inside")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position.inside = c(0.70, 0.25),
      legend.title = element_blank(),
      legend.background = element_blank(),
      title = element_text(size = 10),
    ) +
    scale_color_manual(values = c(
      "HCC Fusion" = "#CC5F5A",
      "HCC Slim" = "#7A848D",
      "AFP" = "#A2AFA6"
    ))
}

roc_plots <- prepared_data %>%
  nest_by(subtype, comparison) %>%
  mutate(
    data = list(data %>% group_by(predictor)),
    plot = list(plot_roc_curve(data, str_c(subtype, ": ", comparison)))
  ) %>%
  select(-data)

roc_p <- reduce(roc_plots$plot, `+`) + plot_layout(nrow = 2)
# tgutil::ggpreview(plot = roc_p, width = 12, height = 6)
ggsave(snakemake@output[[1]], plot = roc_p, width = 12, height = 6)

# PR Curves-----
plot_pr_curve <- function (data, .title) {
  data %>%
    pr_curve(target, value) %>%
    ggplot(aes(recall, precision, color = predictor)) +
    geom_path() +
    coord_equal() +
    ggtitle(.title) +
    lims(x = c(0, 1), y = c(0, 1)) +
    guides(color = guide_legend(position = "inside")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position.inside = c(0.40, 0.25),
      legend.title = element_blank(),
      legend.background = element_blank(),
      title = element_text(size = 10),
    ) +
    scale_color_manual(values = c(
      "HCC Fusion" = "#CC5F5A",
      "HCC Slim" = "#7A848D",
      "AFP" = "#A2AFA6"
    ))
}

pr_plots <- prepared_data %>%
  nest_by(subtype, comparison) %>%
  mutate(
    data = list(data %>% group_by(predictor)),
    plot = list(plot_pr_curve(data, str_c(subtype, ": ", comparison)))
  ) %>%
  select(-data)

pr_p <- reduce(pr_plots$plot, `+`) + plot_layout(nrow = 2)
# tgutil::ggpreview(plot = pr_p, width = 12, height = 6)
ggsave(snakemake@output[[2]], plot = pr_p, width = 12, height = 6)