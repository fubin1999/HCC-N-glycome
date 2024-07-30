# Perform consensus clustering to identify molecular subtypes.
# 1. To avoid the influence of missing values, we will only
# use glycans that are present in at least 80% of the samples.
# 2. Furthermore, to prevent the influence of batch effects, we will
# use the ComBat algorithm to correct for batch effects.

library(tidyverse)
library(sva)
library(ConsensusClusterPlus)

# raw_data <- read_csv("results/data/prepared/raw_abundance.csv")
# processed_data <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# plates <- read_csv("data/plates.csv")

raw_data <- read_csv(snakemake@input[[1]])
processed_data <- read_csv(snakemake@input[[2]])
groups <- read_csv(snakemake@input[[3]])
plates <- read_csv(snakemake@input[[4]])

# 1. Filter out glycans with more than 10% missing values
low_missing_glycans <- raw_data %>%
  left_join(groups, by = "sample") %>%
  filter(group == "HCC") %>%
  select(-group) %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  group_by(glycan) %>%
  summarise(missing = sum(is.na(value)) / n()) %>%
  filter(missing < 0.1) %>%
  pull(glycan)

data_for_cc <- processed_data %>%
  left_join(groups, by = "sample") %>%
  filter(group == "HCC") %>%
  select(-group) %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  filter(glycan %in% low_missing_glycans) %>%
  pivot_wider(names_from = sample, values_from = value) %>%
  column_to_rownames(var = "glycan")
data_for_cc <- log(data_for_cc)
data_for_cc <- as.matrix(data_for_cc)

# 2. Correct for batch effects
plates <- plates %>%
  select(sample, plate) %>%
  filter(sample %in% colnames(data_for_cc))
batches <- plates %>%
  column_to_rownames(var = "sample") %>%
  .[colnames(data_for_cc),]

modcombat <- model.matrix(~1, data = plates)
combat_data <- ComBat(dat = data_for_cc, batch = batches, mod = modcombat)

# 3. Perform consensus clustering
normed <- sweep(combat_data, 1, apply(combat_data, 1, median, na.rm = T))
cc_result <- ConsensusClusterPlus(
  normed, maxK = 8, reps = 1000, pItem = 0.8, pFeature = 1,
  # title = "results/figures/subtypes/cc_result",
  title = snakemake@output[[1]],
  clusterAlg = "hc", distance = "pearson", seed = 123, plot = "png"
)

cluster_classes <- cc_result[[3]][["consensusClass"]]
cluster_classes <- data.frame(sample = names(cluster_classes), class = as.integer(cluster_classes))
write.csv(cluster_classes, snakemake@output[[2]], row.names = FALSE)