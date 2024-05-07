source("renv/activate.R")

library(tidyverse)
library(factoextra)

# Load data-----
#abundance_file <- "results/data/prepared/processed_abundance.csv"
#plates_file <- "data/plates.csv"
abundance_file <- snakemake@input[[1]]
plates_file <- snakemake@input[[2]]
abundance <- read_csv(abundance_file)
plates <- read_csv(plates_file) |> 
  select(sample, plate)

data <- abundance |> 
  left_join(plates, by = "sample")

# PCA on abundance data-----
pca_data <- data |> 
  select(sample, glycan, value) |>
  pivot_wider(names_from = glycan, values_from = value) |> 
  column_to_rownames(var = "sample")
pca <- prcomp(pca_data, center = TRUE, scale = TRUE)

# Plot PCA-----
fviz_pca_ind(
  pca, 
  geom.ind = "point",
  col.ind = data |> distinct(sample, plate) |> pull(plate) |> as.factor(),
  addEllipses = TRUE,
)
# tgutil::ggpreview(width = 6, height = 5)
ggsave(snakemake@output[[1]], width = 6, height = 5)
# ggsave("results/figures/data_quality/batch_effect_pca.pdf", width = 6, height = 5)
file.remove("Rplots.pdf")