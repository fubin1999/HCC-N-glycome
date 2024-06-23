library(tidyverse)

# Read data-----
abundance <- read_csv(snakemake@input[["abundance"]])
groups <- read_csv(snakemake@input[["groups"]])
ancova_result <- read_csv(snakemake@input[["ancova_result"]])

diff_glycans <- ancova_result %>%
  filter(p.adj < 0.05, Effect == "group") %>%
  pull(glycan)

data <- abundance %>%
  inner_join(groups, by = "sample") %>%
  filter(group != "QC") %>%
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) %>%
  pivot_longer(-c(sample, group), names_to = "glycan", values_to = "value") %>%
  mutate(
    log_value = log2(value),
    z_score = (log_value - mean(log_value)) / sd(log_value)
  )

diff_data <- data %>%
  filter(glycan %in% diff_glycans)

# Draw boxplots-----
p <- ggplot(diff_data, aes(group, z_score, fill = group)) +
  geom_violin() +
  facet_wrap(~ glycan, scales = "free", ncol = 6) +
  labs(x = "", y = "z-scores (log abundance)") +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = c("HC" = "#7A848D", "CHB" = "#A2AFA6", "LC" = "#FEC37D", "HCC" = "#CC5F5A"))
# tgutil::ggpreview(width = 12, height = 16)
ggsave(snakemake@output[[1]], plot = p, width = 12, height = 16)