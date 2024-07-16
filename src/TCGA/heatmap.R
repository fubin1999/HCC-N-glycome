library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# load("results/data/TCGA/prepared_data.rda")
# dea_result <- read_csv("results/data/TCGA/dea_results.csv")
# glyco_genes <- read_csv("data/glycogenes.csv")

load(snakemake@input[[1]])
dea_result <- read_csv(snakemake@input[[2]])
glyco_genes <- read_csv(snakemake@input[[3]])

expr_mat <- assay(data, "tpm_unstrand")
gene_info <- rowData(data) %>% as_tibble()
sample_info <- colData(data) %>%
  as_tibble() %>%
  sjmisc::rec(sample_type, rec = "Primary Tumor=PT;Solid Tissue Normal=TN", suffix = "")

gene_info_sub <- gene_info %>%
  semi_join(glyco_genes, by = "gene_name")
expr_mat_sub <- expr_mat[gene_info_sub$gene_id,]
rownames(expr_mat_sub) <- gene_info_sub$gene_name

# Normalize to make TN samples look alike
z_scores <- expr_mat_sub %>%
  as.data.frame() %>%
  rownames_to_column("gene_name") %>%
  as_tibble() %>%
  pivot_longer(-gene_name, names_to = "barcode", values_to = "expr") %>%
  left_join(sample_info %>% select(barcode, sample_type), by = "barcode") %>%
  group_by(gene_name) %>%
  mutate(
    log_expr = log(expr * 1000 + 1),
    z_score = (log_expr - mean(log_expr)) / sd(log_expr)
  ) %>%
  ungroup()

TN_medians <- z_scores %>%
  filter(sample_type == "TN") %>%
  summarise(median_z_score = median(z_score), .by = gene_name)

normed <- z_scores %>%
  left_join(TN_medians, by = "gene_name") %>%
  mutate(normed = z_score - median_z_score) %>%
  select(gene_name, barcode, normed)

z_score_means <- normed %>%
  left_join(sample_info, by = "barcode") %>%
  summarise(mean_z_score = mean(normed), .by = c(gene_name, sample_type))

mat_to_plot <- z_score_means %>%
  pivot_wider(names_from = gene_name, values_from = mean_z_score) %>%
  column_to_rownames("sample_type") %>%
  as.matrix()

col_split <- glyco_genes %>%
  column_to_rownames("gene_name") %>%
  .[gene_info_sub$gene_name, 1]

col_anno <- dea_result %>%
  semi_join(glyco_genes, by = "gene_name") %>%
  column_to_rownames("gene_name") %>%
  .[gene_info_sub$gene_name, c("logFC", "FDR")]

pdf(snakemake@output[[1]], width = 15, height = 2)
col_fun <- colorRamp2(c(-2., 0, 2.), c("#27408B", "white", "#CD0000"))
Heatmap(
  mat_to_plot,
  col = col_fun,
  name = "Z-score",
  border = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_side = "left",
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  column_split = col_split,
  column_gap = unit(2, "mm"),
  cluster_column_slices = FALSE,
  column_names_rot = -60,
  rect_gp = gpar(col = "black", lwd = 1)
)
dev.off()