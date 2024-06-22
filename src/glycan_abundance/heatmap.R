library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Read data-----
abundance <- read_csv("results/data/prepared/processed_abundance.csv")
groups <- read_csv("results/data/prepared/groups.csv")
anova_result <- read_csv("results/data/glycan_abundance/ancova_for_glycans.csv")

abundance <- read_csv(snakemake@input[["abundance"]])
groups <- read_csv(snakemake@input[["groups"]])
anova_result <- read_csv(snakemake@input[["ancova_result"]])

groups <- groups |> 
  filter(group != "QC")
data <- abundance |>
  pivot_longer(-sample, names_to = "glycan", values_to = "value") |>
  semi_join(groups, by = "sample")

# Prepare for heatmap-----
diff_glycans <- anova_result |> 
  filter(p.adj < 0.05, Effect == "group") |>
  pull(glycan)

mat <- data |> 
  filter(glycan %in% diff_glycans) |> 
  mutate(value = log(value)) |> 
  group_by(glycan) |> 
  mutate(value = (value - mean(value)) / sd(value)) |> 
  pivot_wider(names_from = "sample", values_from = "value") |> 
  column_to_rownames("glycan") |> 
  as.matrix()

# Shuffle columns-----
set.seed(42)
mat <- mat[, sample(ncol(mat))]

# Plot heatmap-----
col <- colorRamp2(c(-2, 0, 2), c("#27408B", "white", "#CD0000"))

col_split <- groups |> 
  column_to_rownames("sample")
col_split <- col_split[colnames(mat), ]
col_split <- factor(col_split, levels = c("HC", "CHB", "LC", "HCC"))

pdf(snakemake@output[[1]], width = 5, height = 8)
set.seed(42)
ht <- Heatmap(
  mat,
  name = "Z-score",
  col = col,
  top_annotation = HeatmapAnnotation(group = anno_block(
    gp = gpar(fill = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A")))
  ),
  border = TRUE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  column_split = col_split,
  row_km = 5,
  row_km_repeats = 100,
  row_title = "GCM%s",
  cluster_row_slices = FALSE,
  heatmap_legend_param = list(direction = "horizontal")
)
ht <- draw(ht, heatmap_legend_side = "bottom")
dev.off()

glycan_cluster <- map(row_order(ht), ~ rownames(mat)[.]) |> 
  map(~ tibble(glycan = .x)) |> 
  list_rbind(names_to = "cluster")
write_csv(glycan_cluster, snakemake@output[[2]])
