library(tidyverse)
library(rstatix)

# eigen_glycan <- read_csv("results/data/glycan_coexpr/eigen_glycans.csv")
# subtypes <- read_csv("results/data/subtypes/consensus_cluster_result.csv")

eigen_glycan <- read_csv(snakemake@input[[1]])
subtypes <- read_csv(snakemake@input[[2]])

data <- eigen_glycan %>%
  mutate(gcm = paste0("GCM", cluster), .keep = "unused") %>%
  right_join(subtypes %>% rename(subtype = class), by = "sample") %>%
  mutate(subtype = paste0("Subtype", subtype)) %>%
  relocate(sample, subtype, gcm)

anova_result <- data %>%
  group_by(gcm) %>%
  anova_test(eigen_glycan ~ subtype) %>%
  adjust_pvalue(method = "BH") %>%
  as_tibble() %>%
  select(gcm, F, ges, p.adj)

post_hoc_result <- data %>%
  group_by(gcm) %>%
  tukey_hsd(eigen_glycan ~ subtype) %>%
  as_tibble() %>%
  select(-c(term, null.value))

write_csv(anova_result, snakemake@output[[1]])
write_csv(post_hoc_result, snakemake@output[[2]])