library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Read data-----
# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# anova_result <- read_csv("results/data/diff_analysis/ancova_for_glycans.csv")
# mp_table <- read_csv("results/data/derived_traits/meta_properties.csv")

abundance <- read_csv(snakemake@input[["abundance"]])
groups <- read_csv(snakemake@input[["groups"]])
anova_result <- read_csv(snakemake@input[["ancova_result"]])
mp_table <- read_csv(snakemake@input[["mp_table"]])

groups <- groups |> 
  filter(group != "QC")
data <- abundance |>
  pivot_longer(-sample, names_to = "glycan", values_to = "value") |>
  semi_join(groups, by = "sample")

# Prepare for heatmap-----
diff_glycans <- anova_result |> 
  filter(p.adj < 0.05, Effect == "group") |>
  pull(glycan)

mat <- data |> 
  filter(glycan %in% diff_glycans) |> 
  mutate(value = log(value)) |> 
  group_by(glycan) |> 
  mutate(value = (value - mean(value)) / sd(value)) |> 
  pivot_wider(names_from = "sample", values_from = "value") |> 
  column_to_rownames("glycan") |> 
  as.matrix()

# Shuffle columns-----
set.seed(42)
mat <- mat[, sample(ncol(mat))]

# Plot heatmap-----
col <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

col_split <- groups |> 
  column_to_rownames("sample")
col_split <- col_split[colnames(mat), ]
col_split <- factor(col_split, levels = c("HC", "CHB", "LC", "HCC"))

rela_abund <- data |> 
  group_by(sample) |> 
  mutate(value = value / sum(value) * 100) |> 
  ungroup() |> 
  summarise(median_abund = median(value), .by = glycan) |> 
  mutate(Abund = log(median_abund)) |> 
  select(-median_abund)

row_anno_df <- mp_table |> 
  select(glycan, type, B, nAnt, nF, nS, nG) |> 
  rename(Type = type, Bisect = B) |> 
  sjmisc::rec(Type, rec = "high_mannose=high-mannose;else=copy", suffix = "") |> 
  filter(glycan %in% diff_glycans) |> 
  left_join(rela_abund, by = "glycan") |> 
  column_to_rownames("glycan") |> 
  mutate(
    nAnt = factor(nAnt, ordered = TRUE),
    nF = factor(nF, ordered = TRUE),
    nS = factor(nS, ordered = TRUE),
    nG = factor(nG, ordered = TRUE)
  )
row_anno_df <- row_anno_df[rownames(mat),]
abund_col <- colorRamp2(c(-4, 4), c("white", "pink2"))
row_anno <- rowAnnotation(
  df = row_anno_df,
  annotation_name_side = "top",
  col = list(
    Type = c("complex" = "#FF8C00", "hybrid" = "#FFD700", "high-mannose" = "#228B22"),
    Bisect = c(`TRUE` = "grey30", `FALSE` = "grey95"),
    nAnt = c(`0` = "grey95", `1` = "#CDE6F6", `2` = "#99C6E4", `3` = "#4888B1", `4` = "#0A446A"),
    nF = c(`0` = "grey95", `1` = "#FF8A8A", `2` = "#EF3A3A"),
    nS = c(`0` = "grey95", `1` = "#D6A7E8", `2` = "#B674D0", `3` = "#8131A1", `4` = "#6D158F"),
    nG = c(`0` = "grey95", `1` = "#FFDCA4", `2` = "#FFCC79", `3` = "#DA972A", `4` = "#AC7213"),
    Abund = abund_col
  )
)

pdf(snakemake@output[[1]], width = 7, height = 8)
set.seed(42)
ht <- Heatmap(
  mat,
  name = "Z-score",
  col = col,
  top_annotation = HeatmapAnnotation(group = anno_block(
    gp = gpar(fill = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A")))
  ),
  right_annotation = row_anno,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  column_split = col_split,
  row_km = 4,
  row_km_repeats = 100,
  row_title = "Cluster %s",
  cluster_row_slices = FALSE,
  heatmap_legend_param = list(direction = "horizontal")
)
ht <- draw(ht, heatmap_legend_side = "bottom")
dev.off()

glycan_cluster <- map(row_order(ht), ~ rownames(mat)[.]) |> 
  map(~ tibble(glycan = .x)) |> 
  list_rbind(names_to = "cluster")
write_csv(glycan_cluster, snakemake@output[[2]])
