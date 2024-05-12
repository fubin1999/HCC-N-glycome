library(tidyverse)

# clinical <- read_csv("results/data/prepared/clinical.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
clinical <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

AFP_data <- clinical |> 
  select(sample, AFP) |> 
  inner_join(groups, by = "sample")

calcu_metrics <- function(data, cutoff) {
  data |> 
    mutate(
      AFP_pos = AFP >= {{ cutoff }},
      TP = AFP_pos & group == "C",
      FP = AFP_pos & group != "C",
      TN = !AFP_pos & group != "C",
      FN = !AFP_pos & group == "C"
    ) |> 
    summarize(
      sensitivity = sum(TP) / (sum(TP) + sum(FN)),
      specificity = sum(TN) / (sum(TN) + sum(FP)),
      accuracy = (sum(TP) + sum(TN)) / n(),
    ) |> 
    mutate(youden = sensitivity + specificity - 1)
}

# Combine the results into a tibble
cutoffs <- c(5, 10, 20, 100, 200, 400)
cutoff_df <- map(cutoffs, ~calcu_metrics(AFP_data, .x)) |> 
  bind_rows() |> 
  mutate(cutoff = cutoffs)

write_csv(cutoff_df, snakemake@output[[1]])