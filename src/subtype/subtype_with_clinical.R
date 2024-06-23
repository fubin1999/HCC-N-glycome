library(tidyverse)
library(rstatix)

# clusters <- read_csv("results/data/subtype/cc_result.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")
clusters <- read_csv(snakemake@input[[1]])
clinical <- read_csv(snakemake@input[[2]])

data <- clusters %>%
  left_join(clinical %>% select(-sex, -age), by = "sample") %>%
  mutate(class = factor(class))

wilcox_result <- data %>%
  pivot_longer(-c(sample, class), names_to = "clinical", values_to = "value") %>%
  group_by(clinical) %>%
  wilcox_test(value ~ class) %>%
  adjust_pvalue(method = "BH")

write_csv(wilcox_result, snakemake@output[[1]])