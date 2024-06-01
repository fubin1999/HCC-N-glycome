library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

# load("results/data/TCGA/prepared_data.rda")
load(snakemake@input[[1]])

prepared <- TCGAanalyze_Preprocessing(data)
unlink("PreprocessingOutput.png")
normed <- TCGAanalyze_Normalization(prepared, geneInfo = geneInfoHT)
filtered <- TCGAanalyze_Filtering(normed, method = "quantile", qnt.cut =  0.25)

samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(filtered),
  typesample = c("NT")
)
samplesTP <- TCGAquery_SampleTypes(
  barcode = colnames(filtered), 
  typesample = c("TP")
)

dea_result <- TCGAanalyze_DEA(
  mat1 = filtered[,samplesNT],
  mat2 = filtered[,samplesTP],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)

write_csv(dea_result, snakemake@output[[1]])