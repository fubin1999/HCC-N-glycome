source("renv/activate.R")

library(tidyverse)
library(broom)
library(forestploter)

# train_data <- read_csv("results/data/ml/train_data.csv") %>%
#   select(sample, group, H5N4F1, H4N4S1, H6N5F1S3, H3N4F1)
train_data <- read_csv(snakemake@input[[1]]) %>%
  select(sample, group, H5N4F1, H4N4S1, H6N5F1S3, H3N4F1)

lr <- glm(
  group ~ H5N4F1 + H4N4S1 + H6N5F1S3 + H3N4F1,
  family = binomial, data = train_data
)
lr_summary <- tidy(lr, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  rename(Term = term) %>%
  mutate(`OR (95% CI)` = sprintf(
    "%.2f (%.2f, %.2f)",
    estimate, conf.low, conf.high)
  ) %>%
  mutate(`P-value` = scales::scientific(p.value, digits = 3)) %>%
  mutate(` ` = paste(rep(" ", 20), collapse = " "))

p <- forest(
  lr_summary %>% select(Term, `P-value`, ` `, `OR (95% CI)`),
  est = lr_summary$estimate,
  lower = lr_summary$conf.low,
  upper = lr_summary$conf.high,
  sizes = lr_summary$std.error,
  ci_column = 3,
  ref_line = 1
)
# tgutil::ggpreview(p, width = 6, height = 3)
ggsave(snakemake@output[[1]], plot = p, width = 6, height = 3)