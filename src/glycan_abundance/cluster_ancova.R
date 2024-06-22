library(tidyverse)
library(rstatix)

source("src/utils/emeans_posthoc.R")

# eigen_glycans <- read_csv("results/data/glycan_abundance/eigen_glycans.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")

eigen_glycans <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]])

data <- eigen_glycans %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical %>% select(sample, age, sex), by = "sample") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))

ancova_result <- data %>%
  group_by(cluster) %>%
  anova_test(eigen_glycan ~ age * sex + group, white.adjust = TRUE)

post_hoc_result <- data %>%
  post_hoc(eigen_glycan ~ age * sex + group, group, cluster)

write_csv(ancova_result, snakemake@output[[1]])
write_csv(post_hoc_result, snakemake@output[[2]])