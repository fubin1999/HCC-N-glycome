library(tidyverse)

# Load data-----
abundance <- read_csv(snakemake@input[[1]])
clinical <- read_csv(snakemake@input[[2]])
groups <- read_csv(snakemake@input[[3]])

# Total Abundance Normalization-----
normed <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "abundance") %>%
  mutate(abundance = abundance / sum(abundance, na.rm = TRUE) * 100, .by = sample) %>%
  pivot_wider(names_from = glycan, values_from = abundance)

clinical_selected <- clinical %>%
  select(sample, AST, ALT, GGT, ALB, TBIL, TP, AFP, child_pugh, AAR, ALBI_score)

prepared <- groups %>%
  inner_join(normed, by = "sample") %>%
  inner_join(clinical_selected, by = "sample") %>%
  filter(group != "QC")

# Train-Test Split-----
set.seed(123)
train_data <- prepared %>%
  slice_sample(prop = 0.7, by = group)
test_data <- prepared %>%
  filter(!sample %in% train_data$sample)

# Save Data-----
write_csv(train_data, snakemake@output[[1]])
write_csv(test_data, snakemake@output[[2]])
