library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(cowplot)

# cor_result <- read_csv("results/data/cor_with_clinical/glycan_cor_with_clinical.csv")
# mp_table <- read_csv("results/data/prepared/meta_properties.csv")

cor_result <- read_csv(snakemake@input[[1]])
mp_table <- read_csv(snakemake@input[[2]])

cor_mat <- cor_result %>%
  select(glycan, clinical, cor) %>%
  pivot_wider(names_from = glycan, values_from = cor) %>%
  column_to_rownames("clinical") %>%
  as.matrix()

p_anno_mat <- cor_result %>%
  mutate(signif = if_else(p.adj < 0.05, "*", "")) %>%
  select(glycan, clinical, signif) %>%
  pivot_wider(names_from = glycan, values_from = signif) %>%
  column_to_rownames("clinical") %>%
  as.matrix()

col_anno_df <- mp_table %>%
  mutate(
    glycan = glycan,
    `Complex Type` = type == "complex",
    `Bisecting` = B,
    `High Branching` = nAnt > 2,
    `Fucosylation` = nF > 0,
    `Sialylation` = nS > 0,
    .keep = "none"
  ) %>%
  column_to_rownames("glycan")
col_anno_df[col_anno_df == FALSE] <- "No"
col_anno_df[col_anno_df == TRUE] <- "Yes"
col_anno_df <- col_anno_df[colnames(cor_mat),]
col_anno_col <- c(Yes = "#DDD0C8", No = "grey95")
col_anno <- HeatmapAnnotation(
  df = col_anno_df,
  show_legend = FALSE,
  col = list(
    `Complex Type` = col_anno_col,
    `Bisecting` = col_anno_col,
    `High Branching` = col_anno_col,
    `Fucosylation` = col_anno_col,
    `Sialylation` = col_anno_col
  ),
  simple_anno_size = unit(0.4, "cm")
)

row_split <- tribble(
  ~clinical, ~clinical_type,
  "HBSAG", "Hepatitis\nRelated",
  "HBEAG", "Hepatitis\nRelated",
  "HBEAB", "Hepatitis\nRelated",
  "HBCAB", "Hepatitis\nRelated",
  "HCV", "Hepatitis\nRelated",
  "TP", "Liver\nFunction",
  "TBIL", "Liver\nFunction",
  "AST", "Liver\nFunction",
  "ALT", "Liver\nFunction",
  "AAR", "Liver\nFunction",
  "GGT", "Liver\nFunction",
  "ALB", "Liver\nFunction",
  "AFP", "Cancer\nMarker",
  "CEA", "Cancer\nMarker",
  "CA199", "Cancer\nMarker"
) %>%
  column_to_rownames("clinical") %>%
  .[rownames(cor_mat), "clinical_type"]

col <- colorRamp2(c(-0.6, 0, 0.6), c("#275D87", "white", "#D26F32"))

pdf(snakemake@output[[1]], width = 14, height = 5)
Heatmap(
  cor_mat,
  name = "Spearman's\nrho",
  col = col,
  cluster_rows = FALSE,
  bottom_annotation = col_anno,
  row_split = row_split,
  row_title_rot = 0,
  cell_fun = function(j, i, x, y, w, h, col) {
    grid.text(p_anno_mat[i, j], x, y, hjust = 0.5, vjust = 0.8)
  }
)
dev.off()