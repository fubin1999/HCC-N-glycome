library(tidyverse)
library(broom)
library(forestploter)
library(grid)

# train_data <- read_csv("results/data/ml/train_data.csv")
train_data <- read_csv(snakemake@input[[1]])

lr <- glm(
  group ~ H7N6F1S3 + H5N4 + H3N5 + H4N3S1 + AFP + HBSAG + HBEAG + HBEAB + AST + ALT + GGT + ALB + TBIL + TP,
  family = binomial, data = train_data
)
lr_summary <- tidy(lr, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  dplyr::rename(Factor = term) %>%
  mutate(`OR (95% CI)` = sprintf(
    "%.2f (%.2f, %.2f)",
    estimate, conf.low, conf.high)
  ) %>%
  mutate(`P-value` = scales::scientific(p.value, digits = 3)) %>%
  mutate(` ` = paste(rep(" ", 20), collapse = " "))

p <- forest(
  lr_summary %>% select(Factor, ` `, `P-value`, `OR (95% CI)`),
  est = lr_summary$estimate,
  lower = lr_summary$conf.low,
  upper = lr_summary$conf.high,
  ci_column = 2,
  ref_line = 1
) %>%
  edit_plot(row = 1:4, gp = gpar(fontface = "bold")) %>%
  add_border(part = "header", row = 1, where = "top") %>%
  add_border(part = "header", row = 1, where = "bottom") %>%
  add_border(part = "body", row = 4, where = "bottom")
# tgutil::ggpreview(p, width = 6, height = 5)
ggsave(snakemake@output[[1]], plot = p, width = 6, height = 5)