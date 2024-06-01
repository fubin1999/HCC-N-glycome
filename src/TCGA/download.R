library(tidyverse)
library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
)

query_result <- getResults(query)
write_csv(query_result, snakemake@output[[1]])
GDCdownload(query, directory = snakemake@output[[2]], files.per.chunk = 10)

GDCprepare(
  query,
  directory = "results/data/GDCdata",
  save = TRUE,
  save.filename = snakemake@output[[3]]
)