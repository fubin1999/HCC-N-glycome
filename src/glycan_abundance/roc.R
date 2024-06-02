library(tidyverse)
library(pROC)
library(patchwork)

# Read data-----
# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  filter(group != "QC")

# Perform ROC analysis
calculate_auc <- function (data) {
  data %>%
    group_by(glycan) %>%
    nest() %>%
    mutate(
      roc = map(data, ~ roc(data = .x, response = group, predictor = value, ci = TRUE)),
      auc = map_dbl(roc, ~ .x$auc),
      ci_upper = map_dbl(roc, ~ .x$ci[[3]]),
      ci_lower = map_dbl(roc, ~ .x$ci[[1]])
    ) %>%
    ungroup() %>%
    select(-data)
}

binary_result <- data %>%
  mutate(
    group = if_else(group == "HCC", "HCC", "Control"),
    group = factor(group, levels = c("Control", "HCC"))
  ) %>%
  calculate_auc()

HC_HCC_result <- data %>%
  filter(group %in% c("HC", "HCC")) %>%
  mutate(group = factor(group, levels = c("HC", "HCC"))) %>%
  calculate_auc()

CHB_HCC_result <- data %>%
  filter(group %in% c("CHB", "HCC")) %>%
  mutate(group = factor(group, levels = c("CHB", "HCC"))) %>%
  calculate_auc()

LC_HCC_result <- data %>%
  filter(group %in% c("LC", "HCC")) %>%
  mutate(group = factor(group, levels = c("LC", "HCC"))) %>%
  calculate_auc()

merged_roc_auc <- bind_rows(
  list(
    "Control/HCC" = binary_result,
    "HC/HCC" = HC_HCC_result,
    "CHB/HCC" = CHB_HCC_result,
    "LC/HCC" = LC_HCC_result
  ),
  .id = "comparison"
) %>%
  select(-roc)

write_csv(merged_roc_auc, snakemake@output[[1]])

# Plot ROC curves-----
plot_roc <- function(data, .title) {
  plot_data <- data %>%
    slice_max(auc, n = 10) %>%
    select(glycan, roc) %>%
    deframe() %>%
    ggroc() +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed") +
    labs(x = "Specificity", y = "Sensitivity", color = "Glycan", title = .title) +
    coord_equal() +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )
}

auc_heatmap <- merged_roc_auc %>%
  mutate(comparison = factor(comparison, levels = c("Control/HCC", "HC/HCC", "CHB/HCC", "LC/HCC"))) %>%
  ggplot(aes(reorder(glycan, auc), comparison, fill = auc)) +
  geom_tile() +
  labs(x = "", y = "", fill = "AUC") +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

roc_curves <- plot_roc(binary_result, "Control/HCC") +
  plot_roc(HC_HCC_result, "HC/HCC") +
  plot_roc(CHB_HCC_result, "CHB/HCC") +
  plot_roc(LC_HCC_result, "LC/HCC") +
  plot_layout(ncol = 2)

p <- roc_curves / auc_heatmap + plot_layout(heights = c(6, 1))
# tgutil::ggpreview(width = 10, height = 10)
ggsave(snakemake@output[[2]], plot = p, width = 10, height = 10)