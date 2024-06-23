library(tidyverse)
library(ConsensusClusterPlus)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

HCC_samples <- groups %>%
  filter(group == "HCC") %>%
  pull(sample)

prepared <- abundance %>%
  filter(sample %in% HCC_samples) %>%
  mutate(across(is.numeric, ~ log2(.x + 1))) %>%
  column_to_rownames("sample") %>%
  t()
prepared <- sweep(prepared, 1, apply(prepared, 1, median, na.rm = T))

results <- ConsensusClusterPlus(
  prepared, maxK = 5, reps = 1000, pItem = 0.8, pFeature = 1,
  title = snakemake@output[[1]],
  #title = "results/figures/subtype/cc",
  clusterAlg = "hc", distance = "pearson", seed = 123, plot = "png"
)
cluster_classes <- results[[2]][["consensusClass"]]
cluster_classes <- data.frame(sample = names(cluster_classes), class = as.integer(cluster_classes))
write.csv(cluster_classes, snakemake@output[[2]], row.names = FALSE)