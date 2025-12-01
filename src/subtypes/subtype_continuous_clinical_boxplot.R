library(tidyverse)
library(ggprism)
library(patchwork)

clinical <- read_csv("results/data/prepared/clinical.csv")
subtypes <- read_csv("results/data/subtypes/consensus_cluster_result.csv")
post_hoc_result <- read_csv("results/data/subtypes/subtype_with_continous_clinical_post_hoc_result.csv")

variable_info <- tribble(
  ~variable, ~name, ~full_name, ~log_transform,
  "age", "Age", "Age", FALSE,
  "AST", "AST", "Aspartate Aminotransferase (U/L)", TRUE,
  "ALT", "ALT", "Alanine Aminotransferase (U/L)", TRUE,
  "GGT", "GGT", "Gamma-glutamyltransferase (U/L)", TRUE,
  "ALB", "ALB", "Albumin (g/L)", FALSE,
  "TBIL", "TBIL", "Total Bilirubin (umol/L)", TRUE,
  "TP", "TP", "Total Protein (g/L)", FALSE,
  "AFP", "AFP", "Alpha-fetoprotein (ng/mL)", TRUE,
  "ALBI_score", "ALBI Score", "ALBI Score", FALSE,
  "AAR", "AAR", "De Ritis Ratio" , TRUE,
)

plot_data <- clinical %>%
  select(sample, all_of(variable_info$variable)) %>%
  right_join(subtypes %>% mutate(subtype = paste0("S", class), .keep = "unused"), by = "sample")

y_max <- max(plot_data[["ALB"]])
y_min <- min(plot_data[["ALB"]])

plot_boxplot <- function (data, var) {
  y_max <- max(data[[var]])
  y_min <- min(data[[var]])
  p_data <- post_hoc_result %>%
    filter(variable == var) %>%
    mutate(
      label = if_else(
        p.adj < 0.001,
        scales::label_scientific(digits = 3)(p.adj),
        scales::label_number(accuracy = 0.001)(p.adj)
      )
    )
  if (variable_info %>% filter(variable == var) %>% pull(log_transform)) {
    log_y_max <- log2(y_max)
    log_y_min <- if_else(y_min == 0, 0, log2(y_min))
    p_data <- p_data %>%
      mutate(y.position = 2 ^ ((group2 - group1) * 0.12 * (log_y_max - log_y_min) + log_y_max))
  } else {
    p_data <- p_data %>%
      mutate(y.position = (group2 - group1) * 0.12 * (y_max - y_min) + y_max)
  }

  ggplot(data, aes(subtype, .data[[var]])) +
    geom_boxplot(aes(color = subtype), outlier.alpha = 0, staplewidth = 0.6) +
    geom_jitter(aes(color = subtype), width = 0.2, size = 0.2) +
    add_pvalue(p_data, bracket.shorten = 0.1, bracket.size = 0.4) +
    guides(color = "none") +
    labs(
      x = "Molecular Subtype",
      y = variable_info$full_name[variable_info$variable == var],
      title = variable_info$name[variable_info$variable == var]
    ) +
    scale_color_manual(values = c(S1 = "#1b9e77", S2 = "#d95f02", S3 = "#7570b3")) +
    scale_y_continuous(
      transform = if_else(variable_info %>% filter(variable == var) %>% pull(log_transform), "log2", "identity"),
      expand = expansion(mult = c(0.05, 0.1))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank()
    )
}

plots <- map(variable_info$variable, ~ plot_boxplot(plot_data, .x))
p <- reduce(plots, `+`) + plot_layout(ncol = 5)
single_width <- 2.5
single_height <- 3
final_width <- single_width * 5
final_height <- single_height * (ceiling(nrow(variable_info) / 5))

# tgutil::ggpreview(p, width = final_width, height = final_height)

ggsave(snakemake@output[[1]], p, width = final_width, height = final_height)

write_csv(plot_data, "results/source_data/Supplementary_Figure_15a.csv")
