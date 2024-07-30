library(tidyverse)
library(ComplexHeatmap)

# abundance <- read_csv("results/data/subtypes/batched_corrected.csv")
# clusters <- read_csv("results/data/subtypes/consensus_cluster_result.csv")
# anova_result <- read_csv("results/data/subtypes/anova.csv")
# coexp_modules <- read_csv("results/data/glycan_coexpr/glycan_clusters.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")

abundance <- read_csv(snakemake@input[["abundance"]])
clusters <- read_csv(snakemake@input[["clusters"]])
anova_result <- read_csv(snakemake@input[["anova_result"]])
coexp_modules <- read_csv(snakemake@input[["coexp_modules"]])
clinical <- read_csv(snakemake@input[["clinical"]])

diff_glycans <- anova_result %>%
  filter(p.adj < 0.05) %>%
  pull(glycan)

data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  filter(glycan %in% diff_glycans) %>%
  right_join(clusters, by = "sample") %>%
  group_by(glycan) %>%
  mutate(
    value = log(value),
    value = as.double(scale(value)),
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = "glycan", values_from = "value")

mat <- data %>%
  select(-class) %>%
  column_to_rownames("sample") %>%
  as.matrix() %>%
  t()

selected_clinical <- c("sample", "sex", "age", "AST", "ALT", "GGT", "ALB", "TBIL", "TP", "AFP", "AAR")
col_anno_df <- data %>%
  mutate(subtype = str_glue("Subtype{class}"), subtype = as.factor(subtype)) %>%
  select(sample, subtype) %>%
  left_join(clinical %>% select(all_of(selected_clinical)), by = "sample") %>%
  column_to_rownames("sample") %>%
  .[colnames(mat),]
set.seed(3)
col_anno <- HeatmapAnnotation(df = col_anno_df)

col_split <- as.factor(col_anno_df$subtype)

gcm <- tibble(glycan = diff_glycans) %>%
  left_join(coexp_modules, by = "glycan") %>%
  mutate(
    cluster = replace_na(cluster, 0),
    cluster = str_glue("GCM{cluster}"),
    cluster = as.factor(cluster)
  ) %>%
  column_to_rownames("glycan") %>%
  .[rownames(mat),]
set.seed(3)
row_anno <- rowAnnotation(
  `Glycan\nCo-expression\nModule` = gcm,
  show_annotation_name = FALSE
)

pdf(snakemake@output[[1]], width = 10, height = 10)
Heatmap(
  mat,
  name = "Z-score",
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = col_split,
  row_split = gcm,
  left_annotation = row_anno,
  top_annotation = col_anno,
  cluster_row_slices = FALSE,
)
dev.off()