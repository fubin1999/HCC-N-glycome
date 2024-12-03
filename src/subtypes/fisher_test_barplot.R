library(tidyverse)
library(patchwork)
library(cowplot)

# clinical <- read_csv("results/data/prepared/clinical.csv")
# subtypes <- read_csv("results/data/subtypes/consensus_cluster_result.csv")
# fisher_result <- read_csv("results/data/subtypes/subtype_with_categoric_clinical_fisher_result.csv")

clinical <- read_csv(snakemake@input[["clinical"]])
subtypes <- read_csv(snakemake@input[["subtypes"]])
fisher_result <- read_csv(snakemake@input[["fisher_result"]])

plot_data <- subtypes %>%
  left_join(clinical, by = "sample") %>%
  select(sample, subtype = class, sex, child_pugh, ALBI_stage, TNM_stage) %>%
  mutate(subtype = paste0("S", subtype)) %>%
  mutate(across(-sample, as.factor))

colors <- list(
  sex = c(F = "#FF4081", M = "#39bcdc"),
  child_pugh = c(A = "#ffd299", B = "#FF7F50", C = "#8f1e00"),
  ALBI_stage = c(I = "#ffc7ff", II = "#a364ff", III = "#342a45"),
  TNM_stage = c(I = "#c6ffe6", II = "#61bc84", III = "#2E8B57", IV = "#345e37")
)

p_val_labels <- fisher_result %>%
  mutate(label = if_else(
    p < 0.001,
    paste0("P = ", scales::label_scientific(digits = 3)(p)),
    paste0("P = ", scales::label_number(accuracy = 0.001)(p))
  )) %>%
  select(variable, label) %>%
  deframe()

clinical_names <- c(
  sex = "Sex",
  child_pugh = "Child-Pugh class",
  ALBI_stage = "ALBI stage",
  TNM_stage = "TNM stage"
)

plot_fisher <- function (data, var) {
  data %>%
    group_by(subtype, .data[[var]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    ggplot(aes(subtype, n, fill = .data[[var]])) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_manual(values = colors[[var]]) +
    labs(
      y = "Percentage", fill = clinical_names[[var]],
      title = p_val_labels[[var]]
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 9),
      legend.position = "bottom",
    )
}

plots <- map(names(colors), ~ plot_fisher(plot_data, .x))
combined_plot <- plot_grid(plotlist = plots, nrow = 1)

# tgutil::ggpreview(combined_plot, width = 7.5, height = 2.5)
ggsave(snakemake@output[[1]], combined_plot, width = 7.5, height = 2.5)