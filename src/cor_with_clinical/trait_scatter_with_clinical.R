library(tidyverse)
library(patchwork)

trait_data <- read_csv(snakemake@input[["traits"]])
groups <- read_csv(snakemake@input[["groups"]])
clinical <- read_csv(snakemake@input[["clinical"]]) %>%
  select(-age, -sex)
corr_result_all <- read_csv(snakemake@input[["corr_result_all"]])
corr_result_HCC <- read_csv(snakemake@input[["corr_result_HCC"]])

data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  inner_join(
    pivot_longer(clinical, -sample, names_to = "clinical", values_to = "clinical_value"),
    by = "sample", relationship = "many-to-many"
  )

HCC_data <- data %>%
  filter(group == "HCC")

trait_clinical_pairs_to_plot <- corr_result_all %>%
  filter(abs(cor) > 0.3, p.adj < 0.05)
trait_clinical_pairs_to_plot_HCC <- corr_result_HCC %>%
  filter(abs(cor) > 0.3, p.adj < 0.05)

plot_scatter <- function (data, cor, p, trait, clinical, title) {
  color <- if_else(cor > 0, "#D26F32", "#275D87")
  annotate_x <- if_else(cor > 0, 0.8, 0.1)
  data %>%
    filter(clinical_value > 0) %>%
    ggplot(aes(log2(clinical_value), log2(value)), size = 0.1) +
    geom_point(color = color, alpha = 0.5, shape = 16) +
    geom_smooth(method = "lm", color = color, fill = color, alpha = 0.2) +
    labs(
      x = str_glue("log2({clinical})"),
      y = str_glue("log2({trait})"),
      title = title
    ) +
    theme_classic()
}

plot_data_all <- data %>%
  semi_join(trait_clinical_pairs_to_plot, by = c("trait", "clinical")) %>%
  nest_by(trait, clinical) %>%
  inner_join(trait_clinical_pairs_to_plot, by = c("trait", "clinical"))
plot_data_HCC <- data %>%
  semi_join(trait_clinical_pairs_to_plot_HCC, by = c("trait", "clinical")) %>%
  nest_by(trait, clinical) %>%
  inner_join(trait_clinical_pairs_to_plot_HCC, by = c("trait", "clinical"))
plot_data <- bind_rows(list(All = plot_data_all, HCC = plot_data_HCC), .id = "data_type")

plots <- plot_data %>%
  mutate(plot = list(plot_scatter(data, cor, p.adj, trait, clinical, str_c(data_type, "Data", sep = " ")))) %>%
  select(trait, clinical, data_type, plot)

nrow <- ceiling(sqrt(nrow(plots)))
ncol <- ceiling(nrow(plots) / nrow)
final_p <- reduce(plots$plot, `+`) + plot_layout(nrow = nrow)
# tgutil::ggpreview(plot = final_p, width = ncol * 3, height = nrow * 3)
ggsave(snakemake@output[[1]], plot = final_p, width = ncol * 3, height = nrow * 3)