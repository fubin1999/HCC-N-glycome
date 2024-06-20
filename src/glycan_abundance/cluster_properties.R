library(tidyverse)
library(ComplexHeatmap)

# clusters <- read_csv("results/data/glycan_abundance/glycan_clusters.csv")
# mp_table <- read_csv("results/data/derived_traits/meta_properties.csv")

clusters <- read_csv(snakemake@input[[1]])
mp_table <- read_csv(snakemake@input[[2]])

mp_names <- setdiff(colnames(mp_table), "glycan")
data <- mp_table %>%
  mutate(
    `Complex Type` = type == "complex",
    `Bisecting` = B,
    `High Branching` = nAnt > 2,
    `Fucosylation` = nF > 0,
    `Sialylation` = nS > 0,
  ) %>%
  select(-all_of(mp_names)) %>%
  right_join(clusters, by = "glycan")

mat <- data %>%
  select(-cluster) %>%
  column_to_rownames("glycan") %>%
  t()
mat[mat == FALSE] <- "No"
mat[mat == TRUE] <- "Yes"

col_split <- data$cluster

pdf(snakemake@output[[1]], width = 10, height = 2.8)
Heatmap(
  mat,
  name = "Feature Presence",
  col = c(Yes = "steelblue", No = "grey90"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = col_split,
  column_title = "Cluster %s",
  column_gap = unit(2, "mm"),
  rect_gp = gpar(col = "white", lwd = 2)
)
dev.off()