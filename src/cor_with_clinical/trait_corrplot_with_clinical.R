library(tidyverse)
library(patchwork)

# full_cor_result <- read_csv("results/data/cor_with_clinical/trait_cor_with_clinical.csv")
# sub_cor_result <- read_csv("results/data/cor_with_clinical/trait_cor_with_clinical_per_group.csv")

full_cor_result <- read_csv(snakemake@input[[1]])
sub_cor_result <- read_csv(snakemake@input[[2]])

data <- full_cor_result %>%
  mutate(group = "All") %>%
  bind_rows(sub_cor_result)

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
  "AAR", "Liver Function",
  "GGT", "Liver Function",
  "ALB", "Liver Function",
  "AFP", "Tumor Markers",
  "CEA", "Tumor Markers",
  "CA199", "Tumor Markers"
)


plot_corrplot <- function (cor_result) {
  plot_data <- cor_result %>%
    mutate(
      cor = case_when(
        cor > 0.5 ~ 0.5,
        cor < -0.5 ~ -0.5,
        .default = cor,
      ),
      signif = p.adj < 0.05,
      cor_direction = case_when(
        p.adj < 0.05 & cor > 0 ~ "Positive",
        p.adj < 0.05 & cor < 0 ~ "Negative",
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
      geom_point(aes(size = abs(cor), color = cor_direction), shape = 16, show.legend = TRUE) +
      scale_size_continuous(range = c(0.1, 4), limits = c(0, 0.5), breaks = c(0.1, 0.3, 0.5)) +
      scale_color_manual(
        limits = c("Positive", "Negative", "Not significant"),
        values = c("Positive" = "#D26F32", "Negative" = "#275D87", "Not significant" = "grey")
      ) +
      labs(size = "|Spearman's rho|", color = "Relation") +
      theme_bw()
  }

  function_p <- plot_data %>%
    filter(clinical_type == "Liver Function") %>%
    plot_heatmap("Liver Function") +
    labs(y = "Liver Function") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
    )
  hepatitis_p <- plot_data %>%
    filter(clinical_type == "Hepatitis Related") %>%
    plot_heatmap("Hepatitis Related") +
    labs(y = "Hepatitis") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
    )
  marker_p <- plot_data %>%
    filter(clinical_type == "Tumor Markers") %>%
    plot_heatmap("Tumor Markers") +
    labs(y = "Markers") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  hepatitis_p / function_p / marker_p +
    plot_layout(guides = "collect") &
    coord_fixed() &
    theme(
      legend.position = "top",
      axis.title.x = element_blank(),
      axis.title.y.left = element_text(size = 10),
    ) &
    scale_y_discrete(position = "right")
}

plots <- data %>%
  nest_by(group) %>%
  mutate(plot = list(plot_corrplot(data))) %>%
  select(-data)
# tgutil::ggpreview(plot = plots[["plot"]][[1]], width = 8, height = 6)

all_p <- plots %>% filter(group == "All") %>% pull(plot) %>% .[[1]]
HC_p <- plots %>% filter(group == "HC") %>% pull(plot) %>% .[[1]]
CHB_p <- plots %>% filter(group == "CHB") %>% pull(plot) %>% .[[1]]
LC_p <- plots %>% filter(group == "LC") %>% pull(plot) %>% .[[1]]
HCC_p <- plots %>% filter(group == "HCC") %>% pull(plot) %>% .[[1]]

ggsave(snakemake@output[["all"]], plot = all_p, width = 8, height = 6)
ggsave(snakemake@output[["HC"]], plot = HC_p, width = 8, height = 6)
ggsave(snakemake@output[["CHB"]], plot = CHB_p, width = 8, height = 6)
ggsave(snakemake@output[["LC"]], plot = LC_p, width = 8, height = 6)
ggsave(snakemake@output[["HCC"]], plot = HCC_p, width = 8, height = 6)