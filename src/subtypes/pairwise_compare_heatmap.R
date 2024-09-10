library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# pairwise_result <- read_csv("results/data/subtypes/subtype_with_other_groups_glycan_ttest.csv")
pairwise_result <- read_csv(snakemake@input[[1]])

plot_data <- pairwise_result %>%
  mutate(comparison = paste0(group1, "_" , group2))

fc_mat <- plot_data %>%
  mutate(log_fc = log2(fold_change)) %>%
  pivot_wider(id_cols = glycan, names_from = comparison, values_from = log_fc) %>%
  column_to_rownames("glycan")

p_sig_mat <- plot_data %>%
  mutate(label = if_else(p.adj < 0.05, "*", "")) %>%
  pivot_wider(id_cols = glycan, names_from = comparison, values_from = label) %>%
  column_to_rownames("glycan")

col_anno_df <- tibble(
  `Glycan Subtype` = c(rep("S1", 3), rep("S2", 3), rep("S3", 3)),
  `Compared with` = rep(c("HC", "CHB", "LC"), 3)
)
col_anno <- HeatmapAnnotation(
  df = col_anno_df,
  col = list(
    `Glycan Subtype` = c("S1" = "#1b9e77", "S2" = "#d95f02", "S3" = "#7570b3"),
    `Compared with` = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D")
  ),
  height = unit(8, "mm"),
  simple_anno_size_adjust = TRUE
)

col <- colorRamp2(c(-2, 0, 2), c("#3d5a80", "white", "#ee6c4d"))

ht <- Heatmap(
  fc_mat,
  col = col,
  cluster_columns = FALSE,
  column_split = c(rep("S1", 3), rep("S2", 3), rep("S3", 3)),
  column_gap = unit(2, "mm"),
  column_title = NULL,
  show_column_names = FALSE,
  top_annotation = col_anno,
  rect_gp = gpar(col = "white", lwd = 2),
  width = ncol(fc_mat) * unit(4, "mm") + unit(4, "mm"),
  height = nrow(fc_mat) * unit(4, "mm"),
  cell_fun = function (j, i, x, y, w, h, col) {
    grid.text(p_sig_mat[i, j], x, y - 0.25 * h, gp = gpar(fontsize = 12, col = "black"))
  },
  heatmap_legend_param = list(
    title = "log2(Fold Change)",
    legend_height = unit(40, "mm"),
    title_position = "lefttop-rot"
  )
)
ht <- draw(ht)
w <- ComplexHeatmap:::width(ht)
w <- convertX(w, "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(ht)
h <- convertY(h, "inch", valueOnly = TRUE)

pdf(snakemake@output[[1]], width = w, height = h)
draw(ht)
dev.off()