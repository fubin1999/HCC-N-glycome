library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

# load("results/data/TCGA/prepared_data.rda")
# clusters <- read_csv("results/data/TCGA/consensus_cluster_result.csv")

load(snakemake@input[[1]])
clusters <- read_csv(snakemake@input[[2]])

prepared <- TCGAanalyze_Preprocessing(data)
unlink("PreprocessingOutput.png")
normed <- TCGAanalyze_Normalization(prepared, geneInfo = geneInfoHT)
filtered <- TCGAanalyze_Filtering(normed, method = "quantile", qnt.cut =  0.25)

cluster_1 <- clusters %>%
  filter(class == 1) %>%
  pull(barcode)

cluster_2 <- clusters %>%
  filter(class == 2) %>%
  pull(barcode)

dea_result <- TCGAanalyze_DEA(
  mat1 = filtered[,cluster_1],
  mat2 = filtered[,cluster_2],
  Cond1type = "Cluster1",
  Cond2type = "Cluster2",
  method = "glmLRT"
)
write_csv(dea_result, snakemake@output[[1]])