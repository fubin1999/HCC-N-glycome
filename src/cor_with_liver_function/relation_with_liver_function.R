# This script calculates the relationships between glycans or derived traits
# and liver function-related clinical variables.
# For continuous clinical variables, correlation is performed.
# For categorical clinical variables, ANOVA or t-test were performed.
# All analyses were performed on two levels:
# on the global level using all sample,
# and on the grouped level using samples from each group.

library(tidyverse)
library(rstatix)

# Load and prepare data-----
# glycan_data <- read_csv("results/data/prepared/processed_abundance.csv")
# trait_data <- read_csv("results/data/prepared/filtered_derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")

glycan_data <- read_csv(snakemake@input[[1]])
trait_data <- read_csv(snakemake@input[[2]])
groups <- read_csv(snakemake@input[[3]])
clinical <- read_csv(snakemake@input[[4]])

glycan_long <- glycan_data %>%
  pivot_longer(-sample, names_to = "feature", values_to = "value") %>%
  mutate(value = log(value + 1))

trait_long <- trait_data %>%
  pivot_longer(-sample, names_to = "feature", values_to = "value")

numeric_lf_variables <- c("AST", "ALT", "GGT", "ALB", "TBIL", "TP", "AAR")
categorical_lf_variables <- c("child_pugh", "ALBI_stage")
liver_function_variables <- c(numeric_lf_variables, categorical_lf_variables)

data <- bind_rows(list(glycan = glycan_long, trait = trait_long), .id = "feature_type")

# Correlation Analysis-----
data_with_numeric_clinical <- data %>%
  right_join(
    clinical %>%
      select(sample, all_of(numeric_lf_variables)) %>%
      pivot_longer(-sample, names_to = "clinical_variable", values_to = "clinical_value"),
    by = "sample", relationship = "many-to-many"
  )

global_cor_result <- data_with_numeric_clinical %>%
  group_by(feature_type, feature, clinical_variable) %>%
  cor_test(value, clinical_value, method = "spearman") %>%
  adjust_pvalue(method = "BH") %>%
  select(feature_type, feature, clinical_variable, cor, statistic, p.adj)

grouped_cor_result <- data_with_numeric_clinical %>%
  left_join(groups, by = "sample") %>%
  group_by(group, feature_type, feature, clinical_variable) %>%
  cor_test(value, clinical_value, method = "spearman") %>%
  adjust_pvalue(method = "BH") %>%
  select(group, feature_type, feature, clinical_variable, cor, statistic, p.adj)

# T-test-----
data_with_categoric_clinical <- data %>%
  right_join(
    clinical %>%
      select(sample, all_of(categorical_lf_variables)) %>%
      # Combine Child-Pugh B/C into one category,
      # and combine ALBI II/III into one category.
      mutate(
        child_pugh = if_else(child_pugh == "A", "A", "B/C"),
        ALBI_stage = if_else(ALBI_stage == "I", "I", "II/III")
      ) %>%
      pivot_longer(-sample, names_to = "clinical_variable", values_to = "clinical_value"),
    by = "sample", relationship = "many-to-many"
  )

ttest_with_effsize <- function (data, .formula) {
  ttest_result <- data %>%
    t_test(.formula) %>%
    adjust_pvalue(method = "BH") %>%
    as_tibble() %>%
    select(-c(.y., n1, n2, df, p, p.adj.signif))

  effsize_result <- data %>%
    cohens_d(value ~ clinical_value) %>%
    as_tibble() %>%
    mutate(effsize = -effsize) %>%
    select(-c(.y., n1, n2))

  ttest_result %>%
    left_join(effsize_result, by = c(group_vars(data), "group1", "group2"))
}

global_ttest_result <- data_with_categoric_clinical %>%
  group_by(feature_type, feature, clinical_variable) %>%
  ttest_with_effsize(value ~ clinical_value)

grouped_ttest_result <- data_with_categoric_clinical %>%
  left_join(groups, by = "sample") %>%
  # As HC and CHB don't have enough samples for high level Child-Pugh or ALBI stage,
  # we only perform t-test on LC and HCC samples.
  filter(group %in% c("LC", "HCC")) %>%
  nest_by(group) %>%
  mutate(
    t_test_result = list(
      data %>%
        group_by(feature_type, feature, clinical_variable) %>%
        ttest_with_effsize(value ~ clinical_value)
    )
  ) %>%
  select(-data) %>%
  unnest(t_test_result)

# Save Results-----
write_csv(global_cor_result, snakemake@output[[1]])
write_csv(grouped_cor_result, snakemake@output[[2]])
write_csv(global_ttest_result, snakemake@output[[3]])
write_csv(grouped_ttest_result, snakemake@output[[4]])
