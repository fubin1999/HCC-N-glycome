library(tidyverse)
library(corrplot)
library(rstatix)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")

abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]])


data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  inner_join(groups, "sample") %>%
  inner_join(clinical %>% select(sample, AFP), "sample")

corr_result <- data %>%
  group_by(glycan) %>%
  cor_test(AFP, value, method = "spearman") %>%
  adjust_pvalue()

write_csv(corr_result, snakemake@output[[1]])

cor_mat_1 <- corr_result$cor[1:(nrow(corr_result)/2)] %>% t()
cor_mat_2 <- corr_result$cor[(nrow(corr_result)/2+1):nrow(corr_result)] %>% t()
colnames(cor_mat_1) <- corr_result$glycan[1:(nrow(corr_result)/2)]
colnames(cor_mat_2) <- corr_result$glycan[(nrow(corr_result)/2+1):nrow(corr_result)]
rownames(cor_mat_1) <- "AFP"
rownames(cor_mat_2) <- "AFP"

p_mat_1 <- corr_result$p.adj[1:(nrow(corr_result)/2)] %>% t()
p_mat_2 <- corr_result$p.adj[(nrow(corr_result)/2+1):nrow(corr_result)] %>% t()
colnames(p_mat_1) <- colnames(cor_mat_1)
colnames(p_mat_2) <- colnames(cor_mat_2)
rownames(p_mat_1) <- "AFP"
rownames(p_mat_2) <- "AFP"

plot_cor <- function (cor_mat, p_mat) {
  corrplot(
    cor_mat,
    p.mat = p_mat,
    method = "ellipse",
    tl.col = "black",
    insig = "label_sig",
    pch.cex = 1,
    cl.pos = "n",
  )
}

pdf(snakemake@output[[2]])
plot_cor(cor_mat_1, p_mat_1)
dev.off()

pdf(snakemake@output[[3]])
plot_cor(cor_mat_2, p_mat_2)
dev.off()
