library(tidyverse)
library(rstatix)

# Read data-----
# subtypes <- read_csv("results/data/subtypes/consensus_cluster_result.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")

subtypes <- read_csv(snakemake@input[[1]])
clinical <- read_csv(snakemake@input[[2]])

data <- subtypes %>%
  mutate(subtype = factor(class), .keep = "unused") %>%
  left_join(clinical, by = "sample")

numeric_variables <- c("age", "AST", "ALT", "GGT", "ALB", "TBIL", "TP", "AFP", "ALBI_score", "AAR")
categoric_variables <- c("sex", "child_pugh", "TNM_stage", "ALBI_stage")

# Continuous variables-----
numeric_data <- data %>%
  select(sample, subtype, all_of(numeric_variables)) %>%
  pivot_longer(-c(sample, subtype), names_to = "variable", values_to = "value")

kruskal_result <- numeric_data %>%
  group_by(variable) %>%
  kruskal_test(value ~ subtype) %>%
  adjust_pvalue() %>%
  as_tibble() %>%
  select(-c(.y., n, df, p, method))

post_hoc_result <- numeric_data %>%
  group_by(variable) %>%
  dunn_test(value ~ subtype) %>%
  as_tibble() %>%
  select(-c(.y., n1, n2, p, p.adj.signif))

write_csv(kruskal_result, snakemake@output[[1]])
write_csv(post_hoc_result, snakemake@output[[2]])

# Categoric variables-----
## Global Fisher's exact test-----
categoric_data <- data %>%
  select(sample, subtype, all_of(categoric_variables)) %>%
  mutate(across(all_of(categoric_variables), factor))

subtype_fisher_test <- function (data, var1, var2 = "subtype") {
  data %>%
    group_by(.data[[var1]], .data[[var2]]) %>%
    count() %>%
    pivot_wider(names_from = all_of(var2), values_from = n) %>%
    column_to_rownames(var = var1) %>%
    as.matrix() %>%
    fisher_test(workspace = 2e+08) %>%
    as_tibble() %>%
    select(-p.signif, -n) %>%
    mutate(variable = var1)
}

fisher_result <- map_dfr(categoric_variables, ~ subtype_fisher_test(categoric_data, .x))

write_csv(fisher_result, snakemake@output[[3]])

## Hypergeometric test -----
# Here we used the same principle as the enrichment analysis,
# to test which two subtypes are overlapped.
# This is the post hoc test for the fisher test above.
hypergeometric_test <- function (x, p.adj.method = "fdr") {
  x <- as.matrix(x)
  p_mat <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(ncol(x))) {
      m <- x[i, j]
      n <- sum(x[i, ])
      M <- sum(x[, j])
      N <- sum(x)
      p_mat[i, j] <- phyper(m, M, N - M, n, lower.tail = FALSE)
    }
  }
  p_mat_adj <- matrix(p.adjust(p_mat, p.adj.method), nrow = nrow(x))
  rownames(p_mat_adj) <- rownames(x)
  colnames(p_mat_adj) <- colnames(x)
  p_mat_adj
}

tidy_hypergeometric_test <- function (data, var1, var2) {
  data %>%
    group_by(.data[[var1]], .data[[var2]]) %>%
    count() %>%
    pivot_wider(names_from = all_of(var2), values_from = n) %>%
    column_to_rownames(var = var1) %>%
    as.matrix() %>%
    hypergeometric_test() %>%
    as.data.frame() %>%
    rownames_to_column(var1) %>%
    pivot_longer(-all_of(var1), names_to = var2, values_to = "p_value")
}

hypergeometric_test_result <- list_rbind(map(
  c("child_pugh", "TNM_stage", "ALBI_stage"),
  ~ categoric_data %>%
    tidy_hypergeometric_test(.x, "subtype") %>%
    rename(clinical_subtype = all_of(.x)) %>%
    mutate(clinical_variable = .x, .before = 1)
))

write_csv(hypergeometric_test_result, snakemake@output[[4]])