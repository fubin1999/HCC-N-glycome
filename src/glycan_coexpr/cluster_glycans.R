library(tidyverse)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(circlize)
library(sClust)
library(smotefamily)

# The second version of this method:
# It contains the following steps:
#   1. Rule out variables that are not significantly different between groups.
#   2. Over-sample the data using SMOTE.
#   3. Separate the data into subsets by grouping.
#   4. Perform consensus clustering on each subset.
#   5. Merge the consensus matrices by element-wise minimum.
#   6. Perform spectral clustering on the merged matrix.

# Load data-----
abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

diff_glycans <- read_csv(snakemake@input[[3]]) %>%
  filter(Effect == "group", p.adj < 0.05) %>%
  pull(glycan)


# 0. Prepare data-----
# Remove QC samples, and log-transform the data
data <- abundance %>%
  left_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(
    across(-c(sample, group), log),
    across(-c(sample, group), ~ as.double(scale(.)))
  )

# 1. Rule out variables that are not significantly different between groups-----
data <- data %>%
  select(sample, group, all_of(diff_glycans))

# 2. Over-sample the data using SMOTE-----
X <- data %>%
  select(-group) %>%
  column_to_rownames("sample")
target <- data$group
smote_result <- SMOTE(X, target)$data

# 3. Separate the data into subsets by grouping-----
subsets <- vector("list", length(unique(data$group)))
names(subsets) <- unique(data$group)
for (group in names(subsets)) {
  subsets[[group]] <- smote_result %>%
    filter(class == !!group) %>%
    select(-class) %>%
    as.matrix()
}

# 4. Perform consensus clustering on each subset-----
cc_results <- vector("list", length(subsets))
names(cc_results) <- names(subsets)
dir.create(snakemake@output[[1]], showWarnings = FALSE)
for (group in names(subsets)) {
  cc <- ConsensusClusterPlus(
    subsets[[group]],
    title = file.path(snakemake@output[[1]], group),
    maxK = 8,
    reps = 100,
    pItem = 0.8,
    pFeature = 1,
    plot = "png",
    distance = "euclidean",
    clusterAlg = "km",
    seed = 123
  )
  cc_results[[group]] <- cc
}

n_clusters <- c(HC = 2, CHB = 2, LC = 2, HCC = 3)

# 5. Merge the consensus matrices by element-wise minimum-----
matrices <- vector("list", length(subsets))
names(matrices) <- names(subsets)
for (group in names(subsets)) {
  matrices[[group]] <- cc_results[[group]][[n_clusters[group]]]$consensusMatrix
}
# The final merged matrix is the minimum values for each cell.
# For cell (i, j), the value is the minimum value of (i, j)s in all matrices
merged_matrix <- Reduce(pmin, matrices)
colnames(merged_matrix) <- setdiff(colnames(data), c("sample", "group"))
rownames(merged_matrix) <- colnames(merged_matrix)

# 6. Perform spectral clustering on the merged matrix-----
W <- checking.gram.similarityMatrix(merged_matrix)
eigen_values <- eigen(W)$values
k <- compute.kclust2(eigen_values, 20)
set.seed(123)
sc_result <- VonLuxburgSC(W, K = 5, flagDiagZero = TRUE, verbose = TRUE)
clusters_df <- tibble(
  glycan = rownames(merged_matrix),
  cluster = sc_result$cluster
)

# Save results-----
write_csv(clusters_df, snakemake@output[[2]])
