library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

dea_result <- read_csv("results/data/TCGA/dea_results.csv")

genes <- dea_result %>%
  filter(FDR < 0.01, abs(logFC) > 1) %>%
  pull(gene_name)

ego <- enrichGO(
  gene = genes,
  universe = dea_result$gene_name,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable = TRUE
)
result <- ego@result