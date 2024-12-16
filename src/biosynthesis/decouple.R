library(tidyverse)
library(decoupleR)


# Load data-----
# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# mp_table <- read_csv("results/data/prepared/meta_properties.csv")

abundance <- read_csv(snakemake@input[[1]])
mp_table <- read_csv(snakemake@input[[2]])


# Decouple-----
mat <- abundance %>%
  column_to_rownames("sample") %>%
  as.matrix() %>%
  log2() %>%
  scale() %>%
  t()

networks <- list()

networks[["FUT8"]] <- mp_table %>%
  mutate(mor = if_else(nFc > 0, 1, 0)) %>%
  select(target = glycan, mor)

networks[["MGAT1"]] <- mp_table %>%
  mutate(mor = if_else(type == "high_mannose", 0, 1)) %>%
  select(target = glycan, mor)

networks[["MGAT2"]] <- mp_table %>%
  mutate(mor = if_else(nAnt > 1, 1, 0)) %>%
  select(target = glycan, mor)

networks[["MGAT3"]] <- mp_table %>%
  mutate(mor = as.integer(B)) %>%
  select(target = glycan, mor)

networks[["MGAT4/5"]] <- mp_table %>%
  mutate(mor = case_when(
    nAnt == 3 ~ 1,
    nAnt == 4 ~ 2,
    .default = 0
  )) %>%
  select(target = glycan, mor)

networks[["B4GALT"]] <- mp_table %>%
  mutate(mor = nG) %>%
  select(target = glycan, mor)

networks[["ST3/6GAL"]] <- mp_table %>%
  mutate(mor = nS) %>%
  select(target = glycan, mor)

network <- bind_rows(networks, .id = "source")

res <- decoupleR::run_mlm(mat, network)

res %>%
  select(glycotransferase = source, sample = condition, activity = score) %>%
  pivot_wider(names_from = sample, values_from = activity) %>%
  write_csv(snakemake@output[[1]])
