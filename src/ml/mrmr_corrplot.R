library(tidyverse)
library(corrplot)

abundance <- read_csv("results/data/prepared/processed_abundance.csv")
mrmr_result <- read_csv("results/data/ml/mrmr_result.csv")
abundance <- read_csv(snakemake@input[[1]])
mrmr_result <- read_csv(snakemake@input[[2]])

selected_glycans <- mrmr_result %>%
  slice_head(n = 4) %>%
  pull(glycan)

plot_data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  mutate(value = log(value)) %>%
  filter(glycan %in% selected_glycans) %>%
  pivot_wider(names_from = glycan, values_from = value) %>%
  column_to_rownames("sample")

pdf(snakemake@output[[1]], width = 3, height = 3)
plot_data %>%
  cor() %>%
  corrplot(
    tl.col = "black",
    type = "lower",
    method = "ellipse",
    addCoef.col = "black",
    order = "original",
    addgrid.col = "grey"
  )
dev.off()
