library(tidyverse)
library(rstatix)

# eigen_glycans <- read_csv("results/data/glycan_abundance/eigen_glycans.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv")

eigen_glycans <- read_csv(snakemake@input[[1]])
groups <- read_csv(snakemake@input[[2]])
clinical <- read_csv(snakemake@input[[3]])

cor_result <- eigen_glycans %>%
  inner_join(groups, by = "sample") %>%
  inner_join(
    clinical %>%
      select(-sex, -age) %>%
      pivot_longer(-sample, names_to = "clinical", values_to = "clinical_value"),
    by = "sample", relationship = "many-to-many"
  ) %>%
  group_by(cluster, clinical) %>%
  cor_test(eigen_glycan, clinical_value, method = "spearman") %>%
  adjust_pvalue("p") %>%
  add_significance("p.adj")

plot_data <- cor_result %>%
  select(cluster, clinical, cor, p.adj.signif) %>%
  mutate(
    cluster_name = str_c("GCM", cluster),
    sig_label = if_else(p.adj.signif == "ns", "", p.adj.signif)
  )

p <- ggplot(plot_data, aes(clinical, reorder(cluster_name, desc(cluster)))) +
  geom_tile(aes(fill = cor), color = "white", linewidth = 1) +
  geom_text(aes(label = sig_label)) +
  scale_fill_gradient2(low = "#27408B", mid = "white", high = "#CD0000") +
  labs(y = "Glycan Cluster", size = "|Spearman's rho|", color = "Relation") +
  coord_equal() +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(),
  )
# tgutil::ggpreview(plot = p, width = 6, height = 2.5)

write_csv(cor_result, snakemake@output[[1]])
ggsave(snakemake@output[[2]], plot = p, width = 6.5, height = 2.5)