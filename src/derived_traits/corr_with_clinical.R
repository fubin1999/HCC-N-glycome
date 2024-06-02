library(tidyverse)
library(rstatix)
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

cor_result <- data %>%
  group_by(trait, clinical) %>%
  cor_test(value, clinical_value, method = "spearman") %>%
  as_tibble() %>%
  select(trait, clinical, cor, p)

clinical_type <- tribble(
  ~clinical, ~clinical_type,
  "HBSAG", "Hepatitis Related",
  "HBEAG", "Hepatitis Related",
  "HBEAB", "Hepatitis Related",
  "HBCAB", "Hepatitis Related",
  "HCV", "Hepatitis Related",
  "TP", "Liver Function",
  "TBIL", "Liver Function",
  "AST", "Liver Function",
  "ALT", "Liver Function",
  "GGT", "Liver Function",
  "ALB", "Liver Function",
  "AFP", "Tumor Markers",
  "CEA", "Tumor Markers",
  "CA199", "Tumor Markers"
)

plot_data <- cor_result %>%
  mutate(
    signif = p < 0.05,
    cor_direction = case_when(
      p < 0.05 & cor > 0 ~ "Positive",
      p < 0.05 & cor < 0 ~ "Negative",
      .default = "Not significant"
    ),
    cor_direction = factor(cor_direction, levels = c("Positive", "Negative", "Not significant"))
  ) %>%
  mutate(trait_type = case_when(
    str_ends(trait, "Fc") ~ "Core-Fucosylation",
    str_ends(trait, "Fa") ~ "Arm-Fucosylation",
    str_ends(trait, "S") ~ "Sialylation",
    str_ends(trait, "G") ~ "Galactosylation",
    str_ends(trait, "B") ~ "Bisectioon",
    .default = "Complexity"
  )) %>%
  arrange(trait_type) %>%
  mutate(trait = factor(trait, levels = unique(trait))) %>%
  left_join(clinical_type, by = "clinical")

plot_heatmap <- function (data, y_axis_title) {
  ggplot(data, aes(trait, clinical)) +
    geom_point(aes(size = abs(cor), color = cor_direction), shape = 16) +
    scale_size_continuous(range = c(0.1, 4), limits = c(0, 0.5), breaks = c(0.1, 0.3, 0.5)) +
    scale_color_manual(values = c(Positive = "#D26F32", Negative = "#275D87", `Not significant` = "grey")) +
    labs(size = "|Corr. Coef.|", color = "Relation") +
    coord_fixed(ratio = 1.5) +
    theme_bw()
}

function_p <- plot_data %>%
  filter(clinical_type == "Liver Function") %>%
  plot_heatmap("Liver Function") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
  )
hepatitis_p <- plot_data %>%
  filter(clinical_type == "Hepatitis Related") %>%
  plot_heatmap("Hepatitis Related") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
  )
marker_p <- plot_data %>%
  filter(clinical_type == "Tumor Markers") %>%
  plot_heatmap("Tumor Markers") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p <- hepatitis_p / function_p / marker_p +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  ) &
  scale_y_discrete(position = "right")
# tgutil::ggpreview(width = 8, height = 6)

write_csv(cor_result, snakemake@output[[1]])
ggsave(snakemake@output[[2]], plot = p, width = 8, height = 6)