library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# load("results/data/TCGA/prepared_data.rda")
# clusters <- read_csv("results/data/TCGA/consensus_cluster_result.csv") %>%
#   mutate(class = str_glue("Cluster {class}"))
# glycogenes <- read_csv("data/glycogenes.csv")$gene_name

load(snakemake@input[[1]])
clusters <- read_csv(snakemake@input[[2]]) %>%
  mutate(class = str_glue("Cluster {class}"))
glycogenes <- read_csv(snakemake@input[[3]])$gene_name

expr_mat <- assay(data, "tpm_unstrand")
gene_info <- rowData(data) %>% as_tibble()
sample_info <- colData(data) %>% as_tibble()

genes_to_show <- gene_info %>%
  filter(gene_name %in% glycogenes) %>%
  pull(gene_name)
gene_ids <- gene_info %>%
  filter(gene_name %in% genes_to_show) %>%
  pull(gene_id)
HCC_samples <- sample_info %>%
  filter(shortLetterCode == "TP") %>%
  pull(barcode)

expr_mat_sub <- expr_mat[gene_ids, HCC_samples]
expr_mat_sub <- t(scale(t(log(expr_mat_sub * 1000 + 1))))
rownames(expr_mat_sub) <- column_to_rownames(gene_info, "gene_id")[rownames(expr_mat_sub), "gene_name"]
col_split <- column_to_rownames(clusters, "barcode")[HCC_samples,1]

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("green", "black", "red"))
pdf(snakemake@output[[1]], width = 5, height = 10)
Heatmap(
  expr_mat_sub,
  name = "Z-score",
  col = col_fun,
  border = TRUE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  column_split = col_split,
)
dev.off()