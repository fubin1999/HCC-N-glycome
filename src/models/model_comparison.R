library(tidyverse)
library(easystats)
library(cowplot)


# Load the data-----
# gcm_data <- read_csv("results/data/glycan_coexpr/eigen_glycans.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv", col_select = c(sample, sex, age, ALBI_score))

gcm_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]], col_select = c(sample, sex, age, ALBI_score))

data <- gcm_data %>%
  inner_join(groups, by = "sample") %>%
  inner_join(clinical, by = "sample") %>%
  mutate(
    group = factor(group, levels = c("HC", "CHB", "LC", "HCC")),
    sex = factor(sex, levels = c("M", "F"))
  )


# Model comparison-----
formulas <- list(
  "group" = eigen_glycan ~ sex + age + group,
  "ALBI" = eigen_glycan ~ sex + age + ALBI_score,
  "group + ALBI" = eigen_glycan ~ sex + age + group + ALBI_score,
  "group * ALBI" = eigen_glycan ~ sex + age + group * ALBI_score
)

data_list <- data %>%
  nest_by(cluster) %>%
  mutate(cluster = paste0("GCM", cluster)) %>%
  deframe()

models <- cross2(names(formulas), names(data_list)) %>%
  map_dfr(~ tibble(
    formula_name = .x[[1]],
    gcm = .x[[2]],
    model = list(lm(formulas[[.x[[1]]]], data = data_list[[.x[[2]]]]))
  ))

compare_result <- models %>%
  mutate(model = set_names(model, formula_name)) %>%
  nest_by(gcm) %>%
  mutate(compare = list(compare_performance(data$model))) %>%
  select(-data) %>%
  unnest(compare) %>%
  select(-Model) %>%
  rename(name = Name)

param_result <- models %>%
  mutate(params = map(model, model_parameters, standardize = "smart")) %>%
  select(-model) %>%
  unnest(params)

model_check_plots <- models %>%
  mutate(plot = map(model, check_model)) %>%
  select(-model) %>%
  mutate(plot_name = paste0(gcm, " ~ ", formula_name, ".pdf"))

# Save the results-----
write_csv(compare_result, snakemake@output[[1]])
write_csv(param_result, snakemake@output[[2]])

save_pdf <- function(.p, name) {
  pdf(name, width = 10, height = 10)
  plot(.p)
  dev.off()
}

dir.create(snakemake@output[[3]], showWarnings = FALSE, recursive = TRUE)
walk2(
  model_check_plots$plot, model_check_plots$plot_name,
  ~ save_pdf(.x, file.path(snakemake@output[[3]], .y))
)
