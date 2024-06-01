library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# data_file <- "results/data/ml/model_comparison.csv"
data_file <- snakemake@input[[1]]
data <- read.csv(data_file, row.names = 1) |> 
  arrange(desc(acc_mean))

col_fun <- colorRamp2(c(0.5, 1), colors = c("white", "#D66460"))

mat <- data |> 
  select(acc_mean, auc_mean, f1_mean) |> 
  rename(ACC = acc_mean, AUC = auc_mean, F1 = f1_mean) |>
  as.matrix()

row_anno <- rowAnnotation(
  SD = anno_barplot(
    data$acc_std, 
    width = unit(0.8, "cm"),
    border = FALSE,
    axis_param = list(direction = "reverse")
  ),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 8),
  annotation_name_rot = 0
)

# output_file <- "results/figures/ml/model_comparison_heatmap.pdf"
output_file <- snakemake@output[[1]]
pdf(output_file, width = 4, height = 8)
hm <- Heatmap(
  mat,
  name = "score",
  col = col_fun,
  left_annotation = row_anno,
  border = TRUE,
  rect_gp = gpar(col = "grey40", lwd = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = c("A", rep("B", nrow(data)-1)),
  row_title = NULL,
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 9, col = "grey30"),
  row_names_gp = gpar(fontsize = 8, fontface = c("bold", rep("plain", nrow(data)-1))),
  heatmap_legend_param = list(
    direction = "horizontal", 
    title = "",
    grid_height = unit(2, "mm")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.3f", data[i, j]),
      x, y, 
      gp = gpar(fontsize = 7, col = "grey20", fontface = ifelse(i == 1, "bold", "plain")))
  }
)
draw(hm, heatmap_legend_side = "bottom")
dev.off()
