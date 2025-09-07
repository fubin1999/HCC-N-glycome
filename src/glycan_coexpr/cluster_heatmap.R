library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Read data-----
abundance <- read_csv("results/data/prepared/processed_abundance.csv")
groups <- read_csv("results/data/prepared/groups.csv")
clusters <- read_csv("results/data/glycan_coexpr/glycan_clusters.csv")

abundance <- read_csv(snakemake@input[["abundance"]])
groups <- read_csv(snakemake@input[["groups"]])
clusters <- read_csv(snakemake@input[["clusters"]])

groups <- groups |> 
  filter(group != "QC")
data <- abundance |>
  pivot_longer(-sample, names_to = "glycan", values_to = "value") |>
  semi_join(groups, by = "sample")

# Prepare for heatmap-----
mat <- data |> 
  filter(glycan %in% clusters$glycan) |>
  mutate(value = log(value)) |> 
  group_by(glycan) |> 
  mutate(value = (value - mean(value)) / sd(value)) |> 
  pivot_wider(names_from = "sample", values_from = "value") |> 
  column_to_rownames("glycan") |> 
  as.matrix()

# Shuffle columns-----
set.seed(123)
mat <- mat[, sample(ncol(mat))]

# Plot heatmap-----
col <- colorRamp2(c(-1.5, 0, 1.5), c("#275D87", "white", "#D26F32"))

col_split <- groups |> 
  column_to_rownames("sample")
col_split <- col_split[colnames(mat), ]
col_split <- factor(col_split, levels = c("HC", "CHB", "LC", "HCC"))

row_split <- clusters %>%
  column_to_rownames("glycan")
row_split <- row_split[rownames(mat), ]

pdf(snakemake@output[[1]], width = 4, height = 6.5)
set.seed(42)
ht <- Heatmap(
  mat,
  name = "Z-score",
  col = col,
  top_annotation = HeatmapAnnotation(group = anno_block(
    gp = gpar(
      fill = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A"),
      col = NA
    ))
  ),
  show_column_names = FALSE,
  show_row_names = FALSE,
  cluster_columns = FALSE,
  column_split = col_split,
  row_split = row_split,
  row_title = "GCM%s",
  cluster_row_slices = FALSE,
  heatmap_legend_param = list(direction = "horizontal")
)
ht <- draw(ht, heatmap_legend_side = "bottom")
dev.off()
