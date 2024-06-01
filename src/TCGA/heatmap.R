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
sample_info <- colData(data) %>% as_tibble()

genes_to_show <- dea_result %>%
  filter(gene_name %in% glyco_genes) %>%
  filter(abs(logFC) > 1, FDR < 0.01) %>%
  pull(gene_name)

gene_info_sub <- gene_info %>%
  filter(gene_name %in% genes_to_show)
expr_mat_sub <- expr_mat[gene_info_sub %>% pull(gene_id),]
rownames(expr_mat_sub) <- gene_info_sub %>% pull(gene_name)
expr_mat_sub <- t(scale(t(log(expr_mat_sub * 1000 + 1))))

col_split <- sample_info %>%
  sjmisc::rec(sample_type, rec = "Primary Tumor=PT;Solid Tissue Normal=TN", suffix = "") %>%
  pull(sample_type)

pdf(snakemake@output[[1]], width = 3, height = 4)
Heatmap(
  expr_mat_sub,
  name = "Z-score",
  border = TRUE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_split = col_split,
  top_annotation = HeatmapAnnotation(group = anno_block(
    gp = gpar(fill = c("#7A848D", "HCC" = "#CC5F5A")))
  ),
)
dev.off()