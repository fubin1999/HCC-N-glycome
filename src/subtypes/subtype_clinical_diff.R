library(tidyverse)
library(rstatix)

# clinical <- read_csv("results/data/prepared/clinical.csv")
# clusters <- read_csv("results/data/subtypes/consensus_cluster_result.csv")

clinical <- read_csv(snakemake@input[[1]])
clusters <- read_csv(snakemake@input[[2]])

data <- clinical %>%
  select(-sex) %>%
  right_join(clusters, by = "sample") %>%
  mutate(class = as.factor(class)) %>%
  pivot_longer(cols = -c(sample, class), names_to = "feature", values_to = "value")

kw_result <- data %>%
  group_by(feature) %>%
  kruskal_test(value ~ class) %>%
  adjust_pvalue(method = "BH")

post_hoc_result <- data %>%
  group_by(feature) %>%
  dunn_test(value ~ class)

write_csv(kw_result, snakemake@output[[1]])
write_csv(post_hoc_result, snakemake@output[[2]])