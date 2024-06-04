library(tidyverse)
library(patchwork)

# trait_data <- read_csv("results/data/derived_traits/filtered_derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv") %>%
#   select(-age, -sex)

trait_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]]) %>%
  select(-age, -sex)

data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  inner_join(
    pivot_longer(clinical, -sample, names_to = "clinical", values_to = "clinical_value"),
    by = "sample", relationship = "many-to-many"
  )

HCC_data <- data %>%
  filter(group == "HCC")

plot_scatter <- function (data, .trait, .clinical, .color, .title) {
  data %>%
    filter(trait == .trait, clinical == .clinical) %>%
    filter(clinical_value > 0) %>%
    ggplot(aes(log(clinical_value), log(value))) +
    geom_point(color = .color, alpha = 0.5, shape = 16) +
    geom_smooth(method = "lm", color = .color, fill = .color, alpha = 0.2) +
    labs(x = str_glue("log({.clinical})"), y = str_glue("log({.trait})"), title = .title) +
    theme_classic()
}

plots <- list()
plots[["CFc_AST"]] <- plot_scatter(data, "CFc", "AST", "#D26F32", "All Groups")
plots[["A2SG_AST"]] <- plot_scatter(data, "A2SG", "AST", "#275D87", "All Groups")
plots[["CFa_HCC"]] <- plot_scatter(data, "CFa", "AST", "#D26F32", "All Groups")
plots[["A2SG_AST_HCC"]] <- plot_scatter(HCC_data, "A2SG", "AST", "#275D87", "HCC Samples")
plots[["CFa_AST_HCC"]] <- plot_scatter(HCC_data, "CFa", "AST", "#D26F32", "HCC Samples")

p <- reduce(plots, `+`) + plot_layout(ncol = 3)
# tgutil::ggpreview(plot = p, width = 9, height = 6)
ggsave(snakemake@output[[1]], plot = p, width = 9, height = 6)

