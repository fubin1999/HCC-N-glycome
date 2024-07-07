library(tidyverse)
library(corrplot)

# abundance <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")

abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

QC_samples <- groups %>%
  filter(group == "QC") %>%
  mutate(
    row_num = row_number(),
    new_sample = paste("QC", row_num, sep = "")
  )

data <- abundance %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "value") %>%
  mutate(log_value = log(value)) %>%
  right_join(QC_samples, by = "sample") %>%
  arrange(row_num)

cor_result <- data %>%
  select(sample = new_sample, glycan, log_value) %>%
  pivot_wider(names_from = sample, values_from = log_value) %>%
  column_to_rownames("glycan") %>%
  cor(method = "pearson")

pdf(snakemake@output[[1]], width = 10, height = 10)
corrplot(
  cor_result,
  method = "color",
  type = "lower",
  tl.col = "black",
  tl.srt = 45,
  addgrid.col = 'white',
  col = COL1("Oranges"),
  col.lim = c(0, 1),
  addCoef.col = "white",
  number.cex = 0.7,
)
dev.off()