library(tidyverse)
library(rstatix)

# activities <- read_csv("results/data/biosynthesis/GT_activities.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

activities <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

prepared <- activities %>%
  pivot_longer(-glycotransferase, names_to = "sample", values_to = "activity") %>%
  left_join(groups, by = "sample") %>%
  filter(!group == "QC") %>%
  mutate(group = factor(group, levels = rev(c("HC", "CHB", "LC", "HCC"))))

diff_result <- prepared %>%
  group_by(glycotransferase) %>%
  t_test(activity ~ group, p.adjust.method = "BH") %>%
  as_tibble() %>%
  select(glycotransferase, group1, group2, statistic, p.adj)

eff_size <- prepared %>%
  group_by(glycotransferase) %>%
  cohens_d(activity ~ group) %>%
  as_tibble() %>%
  select(glycotransferase, group1, group2, effsize)

final_result <- diff_result %>%
  left_join(eff_size, by = c("glycotransferase", "group1", "group2"))

write_csv(final_result, snakemake@output[[1]])