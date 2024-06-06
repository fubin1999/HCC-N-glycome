library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# load("results/data/TCGA/prepared_data.rda")
# dea_result <- read_csv("results/data/TCGA/dea_results.csv")
# glyco_genes <- read_csv("data/glycogenes.csv") %>% pull(gene_name)

load(snakemake@input[[1]])
dea_result <- read_csv(snakemake@input[[2]])
glyco_genes <- read_csv(snakemake@input[[3]]) %>% pull(gene_name)

expr_mat <- assay(data, "tpm_unstrand")
gene_info <- rowData(data) %>% as_tibble()
sample_info <- colData(data) %>%
  as_tibble() %>%
  sjmisc::rec(sample_type, rec = "Primary Tumor=PT;Solid Tissue Normal=TN", suffix = "")

genes_to_show <- dea_result %>%
  filter(gene_name %in% glyco_genes) %>%
  pull(gene_name)

gene_info_sub <- gene_info %>%
  filter(gene_name %in% genes_to_show)
expr_mat_sub <- expr_mat[gene_info_sub %>% pull(gene_id),]
rownames(expr_mat_sub) <- gene_info_sub %>% pull(gene_name)

# Normalize to make TN samples look alike
z_scores <- expr_mat_sub %>%
  as.data.frame() %>%
  rownames_to_column("gene_name") %>%
  as_tibble() %>%
  pivot_longer(-gene_name, names_to = "barcode", values_to = "expr") %>%
  group_by(gene_name) %>%
  mutate(
    log_expr = log(expr * 1000 + 1),
    z_score = (log_expr - mean(log_expr)) / sd(log_expr)
  ) %>%
  ungroup()

TN_medians <- z_scores %>%
  left_join(sample_info %>% select(barcode, sample_type), by = "barcode") %>%
  filter(sample_type == "TN") %>%
  summarise(median_z_score = median(z_score), .by = gene_name)

normed <- z_scores %>%
  left_join(TN_medians, by = "gene_name") %>%
  mutate(normed = z_score - median_z_score) %>%
  select(gene_name, barcode, normed)

mat_to_plot <- normed %>%
  pivot_wider(names_from = barcode, values_from = normed) %>%
  column_to_rownames("gene_name") %>%
  as.matrix() %>%
  .[rownames(expr_mat_sub), colnames(expr_mat_sub)]

col_split <- sample_info$sample_type

pdf(snakemake@output[[1]], width = 4, height = 6)
col_fun <- colorRamp2(c(-2., 0, 2.), c("green", "black", "red"))
Heatmap(
  mat_to_plot,
  col = col_fun,
  name = "Z-score",
  border = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  column_split = col_split,
  top_annotation = HeatmapAnnotation(group = anno_block(
    gp = gpar(fill = c("#CC5F5A", "#7A848D")))
  ),
)
dev.off()