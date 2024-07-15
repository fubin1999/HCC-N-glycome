library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# load("results/data/TCGA/prepared_data.rda")
load(snakemake@input[[1]])

sample_info <- colData(data) %>%
  as_tibble()

anno_data <- sample_info %>%
  select(
    barcode,
    `Sample Type` = definition,
    `Weight` = initial_weight,
    `AJCC Pathologic Stage` = ajcc_pathologic_stage,
    `Primary Diagnosis` = primary_diagnosis,
    `Prior Malignancy` = prior_malignancy,
    `Race` = race,
    `Gender` = gender,
    `Ethnicity` = ethnicity,
    `Vital Status` = vital_status
  ) %>%
  arrange(`Sample Type`) %>%
  column_to_rownames("barcode")

set.seed(14)
pdf(snakemake@output[[1]], width = 8, height = 5)
ha <- HeatmapAnnotation(df = anno_data)
ht <- Heatmap(
  matrix(NA, nrow = 0, ncol = nrow(anno_data)),
  top_annotation = ha,
)
draw(ht, annotation_legend_side = "bottom")
dev.off()