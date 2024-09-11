library(tidyverse)
library(rstatix)
library(ComplexHeatmap)
library(circlize)

# Read and prepare data-----
# clinical <- read_csv("results/data/prepared/clinical.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# subtypes <- read_csv("results/data/subtypes/consensus_cluster_result.csv")

clinical <- read_csv(snakemake@input[["clinical"]])
groups <- read_csv(snakemake@input[["groups"]])
subtypes <- read_csv(snakemake@input[["subtypes"]])

numeric_variables <- c("age", "AST", "ALT", "GGT", "ALB", "TBIL", "TP", "AFP", "ALBI_score", "AAR")

prepared <- clinical %>%
  select(sample, all_of(numeric_variables)) %>%
  pivot_longer(-sample, names_to = "clinical_variable", values_to = "value") %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  left_join(subtypes %>% mutate(subtype = paste0("S", class), .keep = "unused"), by = "sample") %>%
  mutate(group = if_else(group == "HCC", subtype, group)) %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "S1", "S2", "S3"))) %>%
  select(-subtype)

# Pairwise Wilcoxon test-----
wilcox_result <- prepared %>%
  group_by(clinical_variable) %>%
  wilcox_test(
    value ~ group,
    comparisons = list(
      c("S1", "HC"), c("S1", "CHB"), c("S1", "LC"),
      c("S2", "HC"), c("S2", "CHB"), c("S2", "LC"),
      c("S3", "HC"), c("S3", "CHB"), c("S3", "LC")
    ),
    p.adjust.method = "none",
  ) %>%
  select(-c(p.adj, p.adj.signif)) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  as_tibble() %>%
  select(-c(.y., n1, n2, p))

effsize_result <- prepared %>%
  group_by(clinical_variable) %>%
  wilcox_effsize(
    value ~ group,
    comparisons = list(
      c("S1", "HC"), c("S1", "CHB"), c("S1", "LC"),
      c("S2", "HC"), c("S2", "CHB"), c("S2", "LC"),
      c("S3", "HC"), c("S3", "CHB"), c("S3", "LC")
    )
  ) %>%
  as_tibble() %>%
  select(-c(.y., n1, n2))

final_result <- wilcox_result %>%
  left_join(effsize_result, by = c("clinical_variable", "group1", "group2"))

write_csv(final_result, snakemake@output[[1]])

# Draw heatmap-----
grouped_median <- prepared %>%
  summarise(median_value = median(value), .by = c("clinical_variable", "group"))

direction <- expand_grid(
  clinical_variable = numeric_variables,
  group1 = c("HC", "CHB", "LC"),
  group2 = c("S1", "S2", "S3")
) %>%
  left_join(grouped_median, by = c("clinical_variable" = "clinical_variable", "group1" = "group")) %>%
  rename(median_value_1 = median_value) %>%
  left_join(grouped_median, by = c("clinical_variable" = "clinical_variable", "group2" = "group")) %>%
  rename(median_value_2 = median_value) %>%
  mutate(direction = if_else(median_value_2 > median_value_1, 1, -1)) %>%
  select(-c(median_value_1, median_value_2))

plot_data <- final_result %>%
  mutate(comparison = paste0(group1, "_" , group2)) %>%
  left_join(direction, by = c("clinical_variable", "group1", "group2")) %>%
  mutate(effsize = effsize * direction)

effsize_mat <- plot_data %>%
  pivot_wider(id_cols = clinical_variable, names_from = comparison, values_from = effsize) %>%
  column_to_rownames("clinical_variable")

p_sig_mat <- plot_data %>%
  mutate(label = if_else(p.adj < 0.05, "*", "")) %>%
  pivot_wider(id_cols = clinical_variable, names_from = comparison, values_from = label) %>%
  column_to_rownames("clinical_variable")

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

col <- colorRamp2(c(-1, 0, 1), c("#3d5a80", "white", "#ee6c4d"))

ht <- Heatmap(
  effsize_mat,
  col = col,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_split = c(rep("S1", 3), rep("S2", 3), rep("S3", 3)),
  column_gap = unit(2, "mm"),
  column_title = NULL,
  show_column_names = FALSE,
  top_annotation = col_anno,
  rect_gp = gpar(col = "white", lwd = 2),
  width = ncol(effsize_mat) * unit(5, "mm") + unit(4, "mm"),
  height = nrow(effsize_mat) * unit(5, "mm"),
  cell_fun = function (j, i, x, y, w, h, col) {
    grid.text(p_sig_mat[i, j], x, y - 0.2 * h, gp = gpar(fontsize = 12, col = "black"))
  },
  heatmap_legend_param = list(
    title = "Effect Size",
    legend_height = unit(40, "mm"),
    title_position = "lefttop-rot"
  )
)
ht <- draw(ht)
w <- ComplexHeatmap:::width(ht)
w <- convertX(w, "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(ht)
h <- convertY(h, "inch", valueOnly = TRUE)

pdf(snakemake@output[[2]], width = w, height = h)
draw(ht)
dev.off()