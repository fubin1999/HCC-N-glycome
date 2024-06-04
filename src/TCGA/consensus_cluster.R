library(SummarizedExperiment)
library(ConsensusClusterPlus)

# load("results/data/TCGA/prepared_data.rda")
# glycogenes <- read.csv("data/glycogenes.csv")$gene_name

load(snakemake@input[[1]])
glycogenes <- read.csv(snakemake@input[[2]])$gene_name

expr_mat <- assay(data, "tpm_unstrand")
sample_info <- colData(data)
gene_info <- rowData(data)
TP_samples <- sample_info[sample_info$shortLetterCode == "TP",]$barcode
glycogene_ids <- gene_info[gene_info$gene_name %in% glycogenes,]$gene_id
sub_expr_mat <- expr_mat[glycogene_ids,TP_samples]
sub_expr_mat <- sweep(sub_expr_mat, 1, apply(sub_expr_mat, 1, median, na.rm = T))

results <- ConsensusClusterPlus(
  sub_expr_mat, maxK = 5, reps = 1000, pItem = 0.8, pFeature = 1,
  title = snakemake@output[[1]],
  # title = "results/figures/TCGA/consensus_cluster/",
  clusterAlg = "hc", distance = "pearson", seed = 123, plot = "png"
)

cluster_classes <- results[[2]][["consensusClass"]]
cluster_classes <- data.frame(barcode = names(cluster_classes), class = as.integer(cluster_classes))
write.csv(cluster_classes, snakemake@output[[2]], row.names = FALSE)