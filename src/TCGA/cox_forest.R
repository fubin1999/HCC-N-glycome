library(tidyverse)
library(forestploter)

# cox_result <- read_csv("results/data/TCGA/cox.csv")
cox_result <- read_csv(snakemake@input[[1]])

prepared <- cox_result %>%
  filter(p.value.adj < 0.05, abs(log(estimate)) > log(1.01)) %>%
  filter(gene_name != "ST3GAL3") %>%
  mutate(`HR (95% CI)` = sprintf(
    "%.2f (%.2f, %.2f)",
    estimate, conf.low, conf.high)
  ) %>%
  mutate(` ` = paste(rep(" ", 20), collapse = " ")) %>%
  mutate(p.value.adj = scales::scientific(p.value.adj, digits = 3)) %>%
  select(
    Gene = gene_name, ` `, `HR (95% CI)`, `Q-value` = p.value.adj,
    estimate, conf.low, conf.high
  )

p <- forest(
  prepared %>% select(-c(estimate, conf.low, conf.high)),
  est = prepared$estimate,
  lower = prepared$conf.low,
  upper = prepared$conf.high,
  ci_column = 2,
  ref_line = 1,
  xlim = c(0.8, 1.2),
  ticks_at = c(0.8, 1, 1.2)
)
# tgutil::ggpreview(plot = p, width = 6, height = 6)
ggsave(snakemake@output[[1]], plot = p, width = 6, height = 6)