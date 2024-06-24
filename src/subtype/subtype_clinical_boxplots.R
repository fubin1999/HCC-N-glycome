library(tidyverse)
library(ggprism)
library(ggsci)
library(patchwork)


# clinical <- read_csv("results/data/prepared/clinical.csv")
# subtypes <- read_csv("results/data/subtype/cc_result.csv")
# diff_result <- read_csv("results/data/subtype/subtype_clinical_diff.csv")

clinical <- read_csv(snakemake@input[[1]])
subtypes <- read_csv(snakemake@input[[2]])
diff_result <- read_csv(snakemake@input[[3]])

diff_clinical <- diff_result %>%
  filter(p.adj < 0.05) %>%
  pull(clinical)

data <- subtypes %>%
  mutate(
    class = str_glue("HCC_S{class}"),
    class = factor(class, levels = c("HCC_S1", "HCC_S2"))
  ) %>%
  left_join(clinical, by = "sample") %>%
  select(sample, class, all_of(diff_clinical)) %>%
  pivot_longer(-c(sample, class), names_to = "clinical", values_to = "value")

y_pos_data <- data %>%
  filter(value > 0) %>%
  mutate(value = log(value)) %>%
  summarise(
    range = max(value) - min(value),
    max = max(value),
    .by = clinical
  )

p_data <- diff_result %>%
  filter(clinical %in% diff_clinical) %>%
  left_join(y_pos_data, by = "clinical") %>%
  mutate(
    label = str_c("p = ", scales::scientific(p.adj)),
    y.position = max + range * 0.05
  ) %>%
  select(clinical, group1, group2, label, y.position)

plot_boxplots <- function (value_data, p_data, .name) {
  value_data %>%
    filter(value > 0) %>%
    ggplot(aes(class, log(value))) +
    geom_boxplot(aes(color = class), outlier.alpha = 0) +
    geom_jitter(aes(color = class), width = 0.25, size = 1) +
    scale_color_npg() +
    add_pvalue(
      data = p_data,
      bracket.size = 0.5,
      label.size = 3.5,
    ) +
    guides(color = "none") +
    labs(y = str_glue("log({.name})")) +
    theme_classic() +
    theme(
      axis.title.x = element_blank()
    )
}

plots <- data %>%
  nest_by(clinical, .key = "value_data") %>%
  left_join(p_data %>% nest_by(clinical, .key = "p_data"), by = "clinical") %>%
  mutate(plot = list(plot_boxplots(value_data, p_data, clinical))) %>%
  select(clinical, plot)

final_p <- reduce(plots$plot, `+`)
# tgutil::ggpreview(width = 7.5, height = 4)
ggsave(snakemake@output[[1]], plot = final_p, width = nrow(plots) * 2.5, height = 4)