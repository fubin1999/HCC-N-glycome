source("renv/activate.R")

library(tidyverse)
library(pROC)

# Read data-----
# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  filter(group != "QC")

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
    select(-data, -roc)
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

merged_result <- bind_rows(
  list(
    "Control/HCC" = binary_result,
    "HC/HCC" = HC_HCC_result,
    "CHB/HCC" = CHB_HCC_result,
    "LC/HCC" = LC_HCC_result
  ),
  .id = "comparison"
)

write_csv(merged_result, snakemake@output[[1]])