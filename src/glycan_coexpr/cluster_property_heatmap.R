library(tidyverse)
library(ComplexHeatmap)

# clusters <- read_csv("results/data/glycan_coexpr/glycan_clusters.csv")
# mp_table <- read_csv("results/data/prepared/meta_properties.csv")

clusters <- read_csv(snakemake@input[[1]])
mp_table <- read_csv(snakemake@input[[2]])

mp_names <- setdiff(colnames(mp_table), "glycan")
data <- mp_table %>%
  mutate(
    `Complex Type` = type == "complex",
    `Bisecting` = B,
    `High Branching` = nAnt > 2,
    `Fucosylation` = nF > 0,
    `Sialylation` = nS > 0,
  ) %>%
  select(-all_of(mp_names)) %>%
  right_join(clusters, by = "glycan") %>%
  # Convert glycan names to 4 numbers: H, N, F, S
  mutate(
    H = as.integer(str_extract(glycan, "H(\\d+)", group = 1)),
    N = as.integer(str_extract(glycan, "N(\\d+)", group = 1)),
    F = as.integer(str_extract(glycan, "F(\\d+)", group = 1)),
    S = as.integer(str_extract(glycan, "S(\\d+)", group = 1))
  ) %>%
  mutate(across(c(H, N, F, S), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(glycan = paste0(H, N, F, S)) %>%
  select(-c(H, N, F, S))

mat <- data %>%
  select(-cluster) %>%
  column_to_rownames("glycan") %>%
  t()
mat[mat == FALSE] <- "No"
mat[mat == TRUE] <- "Yes"

col_split <- data$cluster

ht <- Heatmap(
  mat,
  name = "Feature Presence",
  col = c(Yes = "steelblue", No = "grey90"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = col_split,
  column_title = "GCM%s",
  column_gap = unit(2, "mm"),
  rect_gp = gpar(col = "white", lwd = 2),
  width = ncol(mat) * unit(4, "mm") + unit(8, "mm"),
  height = nrow(mat) * unit(4, "mm")
)
ht <- draw(ht)
w <- ComplexHeatmap:::width(ht)
w <- convertX(w, "inch", valueOnly = TRUE)
h <- ComplexHeatmap:::height(ht)
h <- convertY(h, "inch", valueOnly = TRUE)

pdf(snakemake@output[[1]], width = w, height = h)
draw(ht)
dev.off()