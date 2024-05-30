source("renv/activate.R")

library(tidyverse)
library(rstatix)
library(ggsignif)

# trait_data <- read_csv("results/data/derived_traits/filtered_derived_traits.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")

trait_data <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]])

data <- trait_data %>%
  pivot_longer(-sample, names_to = "trait", values_to = "value") %>%
  inner_join(groups, by = "sample") %>%
  filter(group == "HCC") %>%
  inner_join(clinical %>% select(sample, AFP, age, sex), by = "sample") %>%
  mutate(AFP_pos = AFP > 10)

ancova_result <- data %>%
  group_by(trait) %>%
  anova_test(value ~ AFP_pos + age + sex, white.adjust = TRUE) %>%
  adjust_pvalue() %>%
  as_tibble()

p <- ancova_result %>%
  filter(trait == "CFc", Effect == "AFP_pos") %>%
  pull(p.adj)

CFc_boxplot <- data %>%
  filter(trait == "CFc") %>%
  mutate(x = if_else(AFP_pos, "AFP Pos.", "AFP Neg.")) %>%
  ggplot(aes(x, value)) +
  geom_boxplot(aes(color = x), outlier.alpha = 0) +
  geom_jitter(aes(color = x), shape = 16, alpha = 0.5) +
  geom_signif(comparisons = list(c("AFP Pos.", "AFP Neg.")),
              annotations = scales::number(p, accuracy = 0.001)) +
  guides(color = "none") +
  labs(x = "", y = "Proportion of core-fucosylated") +
  theme_classic() +
  scale_color_manual(values = c("#F49893", "#A63731")) +
  scale_y_continuous(labels = scales::percent_format())
# tgutil::ggpreview(width = 3, height = 4)

write_csv(ancova_result, snakemake@output[[1]])
ggsave(snakemake@output[[2]], plot = CFc_boxplot, width = 3, height = 4)
