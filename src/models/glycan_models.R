library(tidyverse)
library(easystats)

# Load the data-----
# data <- read_csv("results/data/prepared/processed_abundance.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")

data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]])

prepared <- data %>%
  pivot_longer(-sample, names_to = "glycan", values_to = "abundance") %>%
  mutate(log_abundance = log(abundance)) %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical %>% select(sample, sex, age, ALBI_score), by = "sample") %>%
  mutate(
    sex = factor(sex, levels = c("M", "F")),
    group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))
  )

# Fit the models-----
models <- prepared %>%
  nest_by(glycan) %>%
  mutate(model = list(lm(log_abundance ~ sex + age + ALBI_score + group, data = data))) %>%
  select(-data)

# Get model performance-----
model_performances <- models %>%
  mutate(performance = list(performance(model))) %>%
  select(-model) %>%
  unnest(performance)

# Get model coefficients-----
model_coefficients <- models %>%
  mutate(coefficients = list(model_parameters(model, standardize = "smart"))) %>%
  select(-model) %>%
  unnest(coefficients) %>%
  mutate(Parameter = case_match(
    Parameter,
    "sexF" ~ "Sex: F",
    "groupLC" ~ "Group: LC",
    "groupHCC" ~ "Group: HCC",
    "groupCHB" ~ "Group: CHB",
    "ALBI_score" ~ "ALBI Score",
    "age" ~ "Age",
    "(Intercept)" ~ "Intercept"
  ))

# Save the results-----
write_csv(model_performances, snakemake@output[[1]])
write_csv(model_coefficients, snakemake@output[[2]])
