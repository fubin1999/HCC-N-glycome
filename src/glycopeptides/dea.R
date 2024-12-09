library(tidyverse)
library(limma)

# data <- read_csv("results/data/glycoproteomics/prepared.csv")
data <- read_csv(snakemake@input[[1]])

expr_mat <- data %>%
  mutate(glycopeptide = paste(peptide, gly_site, glycan_composition, sep = "_")) %>%
  select(glycopeptide, sample, log2_ratio) %>%
  pivot_wider(names_from = sample, values_from = log2_ratio) %>%
  column_to_rownames(var = "glycopeptide")

groups <- factor(
  rep(c("HC", "CHB", "LC", "HCC"), each = 4),
  levels = c("HC", "CHB", "LC", "HCC")
)

design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
rownames(design) <- colnames(expr_mat)

contrasts <- c("HC-CHB", "HC-LC", "HC-HCC")

limma_dea <- function(expr_mat, design, contrast) {
  contrast.matrix <- makeContrasts(contrasts = contrast, levels = design)
  fit <- lmFit(expr_mat, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  topTable(fit2, coef = 1, n = Inf) %>%
    rownames_to_column(var = "glycopeptide") %>%
    as_tibble() %>%
    mutate(contrast = contrast)
}

results <- map(contrasts, ~ limma_dea(expr_mat, design, .x)) %>%
  list_rbind() %>%
  separate_wider_delim(contrast, "-", names = c("group1", "group2")) %>%
  relocate(glycopeptide, group1, group2, everything())

write_csv(results, snakemake@output[[1]])