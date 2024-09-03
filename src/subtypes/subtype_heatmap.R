library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
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

col_anno_df <- data %>%
  mutate(subtype = str_glue("Subtype {class}"), subtype = as.factor(subtype)) %>%
  select(sample, subtype) %>%
  left_join(clinical %>% select(sample, child_pugh, ALBI_stage, TNM_stage), by = "sample") %>%
  rename(
    `Glycan Subtype` = subtype,
    `Child-Pugh Class` = child_pugh,
    `ALBI Stage` = ALBI_stage,
    `TNM Stage` = TNM_stage,
  ) %>%
  column_to_rownames("sample") %>%
  .[colnames(mat),]
set.seed(2)
col_anno <- HeatmapAnnotation(
  df = col_anno_df,
  col = list(`Glycan Subtype` = c(`Subtype 1` = "#1b9e77", `Subtype 2` = "#d95f02", `Subtype 3` = "#7570b3")),
  annotation_height = unit(0.5, "cm")
)

col_split <- as.factor(col_anno_df[["Glycan Subtype"]])

gcm <- tibble(glycan = diff_glycans) %>%
  left_join(coexp_modules, by = "glycan") %>%
  mutate(
    cluster = replace_na(cluster, 0),
    cluster = str_glue("GCM{cluster}"),
    cluster = as.factor(cluster)
  ) %>%
  column_to_rownames("glycan") %>%
  .[rownames(mat),]
set.seed(15)
row_anno <- rowAnnotation(
  `Glycan\nCo-expression\nModule` = gcm,
  show_annotation_name = FALSE
)

col <- colorRamp2(c(-2, 0, 2), c("#275D87", "white", "#D26F32"))

pdf(snakemake@output[[1]], width = 6, height = 5)
set.seed(1)
Heatmap(
  mat,
  name = "Z-score",
  col = col,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = col_split,
  row_km = 4,
  row_km_repeats = 100,
  row_title = NULL,
  left_annotation = row_anno,
  top_annotation = col_anno,
  cluster_row_slices = FALSE
)
dev.off()