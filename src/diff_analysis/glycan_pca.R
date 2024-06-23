library(tidyverse)
library(broom)
library(patchwork)

abundance <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])

data <- groups %>%
  inner_join(abundance, by = "sample") %>%
  filter(group != "QC")

pca_fit <- data %>%
  select(-sample, -group) %>%
  prcomp(scale = TRUE)

plot_data <- pca_fit %>%
  augment(data) %>%
  select(sample, group, PC1 = .fittedPC1, PC2 = .fittedPC2)

HC_p <- plot_data %>%
  sjmisc::rec(group, rec = "CHB,LC=Others;else=copy", suffix = "") %>%
  mutate(group = factor(group, levels = c("HC", "HCC", "Others"))) %>%
  ggplot(aes(PC2, PC1, color = group, alpha = group)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c("HC" = "#7A848D", "HCC" = "#CC5F5A", "Others" = "black")) +
  scale_alpha_manual(values = c("HC" = 1, "HCC" = 1, "Others" = 0.1)) +
  theme_classic()

MC_p <- plot_data %>%
  sjmisc::rec(group, rec = "HC,LC=Others;else=copy", suffix = "") %>%
  mutate(group = factor(group, levels = c("CHB", "HCC", "Others"))) %>%
  ggplot(aes(PC2, PC1, color = group, alpha = group)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c("CHB" = "#A2AFA6", "HCC" = "#CC5F5A", "Others" = "black")) +
  scale_alpha_manual(values = c("CHB" = 1, "HCC" = 1, "Others" = 0.1)) +
  theme_classic()

YC_p <- plot_data %>%
  sjmisc::rec(group, rec = "HC,CHB=Others;else=copy", suffix = "") %>%
  mutate(group = factor(group, levels = c("LC", "HCC", "Others"))) %>%
  ggplot(aes(PC2, PC1, color = group, alpha = group)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c("LC" = "#FEC37D", "HCC" = "#CC5F5A", "Others" = "black")) +
  scale_alpha_manual(values = c("LC" = 1, "HCC" = 1, "Others" = 0.1)) +
  theme_classic()

p <- HC_p + MC_p + YC_p & theme(
  legend.position = "bottom"
)
# tgutil::ggpreview(width = 10, height = 4)
ggsave(snakemake@output[[1]], plot = p, width = 10, height = 4)