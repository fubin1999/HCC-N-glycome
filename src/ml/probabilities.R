library(tidyverse)


preds <- read_csv("results/data/ml/preds.csv")
clinical <- read_csv("results/data/prepared/unfiltered_clinical.csv")
groups <- read_csv("results/data/prepared/unfiltered_groups.csv")


prepared <- clinical %>% 
  select(sample, AFP, child_pugh, ALBI_stage, TNM_stage) %>% 
  inner_join(
    preds %>% 
      filter(dataset == "test", model == "global") %>% 
      select(-dataset, -model), 
    by = "sample"
  ) %>% 
  left_join(groups, by = "sample") %>% 
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))


plots <- list()


plots[["group"]] <- ggplot(prepared, aes(group, proba)) +
  ggbeeswarm::geom_quasirandom(aes(color = group), width = 0.4) +
  ggsignif::geom_signif(
    comparisons = list(
      c("HC", "HCC"), 
      c("CHB", "HCC"), 
      c("LC", "HCC"), 
      c("HC", "LC"), 
      c("CHB", "LC")
    ),
    step_increase = 0.17, vjust = -0.2, textsize = 3, tip_length = 0,
    map_signif_level = function(p) {
      sci_num <- scales::label_scientific(digits = 3)
      paste0("p = ", sci_num(p))
    }
  ) +
  scale_color_manual(values = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A")) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.1)),
  ) +
  labs(x = "Group", y = "Predicted Probability") +
  guides(color = "none") +
  theme_classic()


plots[["ALBI"]] <- ggplot(prepared, aes(ALBI_stage, proba)) +
  ggbeeswarm::geom_quasirandom(aes(color = ALBI_stage), width = 0.4) +
  ggsignif::geom_signif(
    comparisons = list(c("I", "II"), c("I", "III")),
    step_increase = 0.13, vjust = -0.2, textsize = 3, tip_length = 0,
    map_signif_level = function(p) {
      sci_num <- scales::label_scientific(digits = 3)
      paste0("p = ", sci_num(p))
    }
  ) +
  scale_color_manual(values = c("#BDD7E7", "#6BAED6", "#08519C")) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.1)),
  ) +
  labs(x = "ALBI Stage", y = "Predicted Probability") +
  guides(color = "none") +
  theme_classic()


plots[["ALBI_HCC"]] <- prepared %>% 
  filter(group == "HCC") %>% 
  ggplot(aes(ALBI_stage, proba)) +
  ggbeeswarm::geom_quasirandom(aes(color = ALBI_stage), width = 0.4) +
  ggsignif::geom_signif(
    comparisons = list(c("I", "II"), c("I", "III")),
    step_increase = 0.13, vjust = -0.2, textsize = 3, tip_length = 0,
    map_signif_level = function(p) {
      sci_num <- scales::label_scientific(digits = 3)
      paste0("p = ", sci_num(p))
    }
  ) +
  scale_color_manual(values = c("#BDD7E7", "#6BAED6", "#08519C")) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.1)),
  ) +
  labs(x = "ALBI Stage (HCC)", y = "Predicted Probability") +
  guides(color = "none") +
  theme_classic()


plots[["child_pugh"]] <- prepared %>% 
  mutate(child_pugh = if_else(child_pugh == "A", "A", "B/C")) %>%
  ggplot(aes(child_pugh, proba)) +
  ggbeeswarm::geom_quasirandom(aes(color = child_pugh), width = 0.4) +
  ggsignif::geom_signif(
    comparisons = list(c("A", "B/C")),
    step_increase = 0.13, vjust = -0.2, textsize = 3, tip_length = 0,
    map_signif_level = function(p) {
      sci_num <- scales::label_scientific(digits = 3)
      paste0("p = ", sci_num(p))
    }
  ) +
  scale_color_manual(values = c("#BDD7E7", "#08519C")) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.1)),
  ) +
  labs(x = "Child-Pugh", y = "Predicted Probability") +
  guides(color = "none") +
  theme_classic()


plots[["child_pugh_HCC"]] <- prepared %>% 
  filter(group == "HCC") %>%
  mutate(child_pugh = if_else(child_pugh == "A", "A", "B/C")) %>%
  ggplot(aes(child_pugh, proba)) +
  ggbeeswarm::geom_quasirandom(aes(color = child_pugh), width = 0.4) +
  ggsignif::geom_signif(
    comparisons = list(c("A", "B/C")),
    step_increase = 0.13, vjust = -0.2, textsize = 3, tip_length = 0,
    map_signif_level = function(p) {
      sci_num <- scales::label_scientific(digits = 3)
      paste0("p = ", sci_num(p))
    }
  ) +
  scale_color_manual(values = c("#BDD7E7", "#08519C")) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.1)),
  ) +
  labs(x = "Child-Pugh (HCC)", y = "Predicted Probability") +
  guides(color = "none") +
  theme_classic()


plots[["TNM"]] <- prepared %>% 
  filter(group == "HCC") %>% 
  ggplot(aes(TNM_stage, proba)) +
  ggbeeswarm::geom_quasirandom(aes(color = TNM_stage), width = 0.4) +
  ggsignif::geom_signif(
    comparisons = list(c("I", "III"), c("I", "IV"), c("II", "IV")),
    step_increase = 0.13, vjust = -0.2, textsize = 3, tip_length = 0,
    map_signif_level = function(p) {
      sci_num <- scales::label_scientific(digits = 3)
      paste0("p = ", sci_num(p))
    }
  ) +
  scale_color_manual(values = RColorBrewer::brewer.pal(5, "Blues")[2:5]) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.1)),
  ) +
  labs(x = "TNM Stage", y = "Predicted Probability") +
  guides(color = "none") +
  theme_classic()


all_p <- cowplot::plot_grid(plotlist = plots, nrow = 1)
tgutil::ggpreview(all_p, width = 12, height = 2.5)
ggsave("results/figures/ml/probabilities.pdf", all_p, width = 12, height = 2.5)
