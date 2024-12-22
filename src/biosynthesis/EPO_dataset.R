library(tidyverse)
library(glymotif)
library(decoupleR)


# data <- readxl::read_excel("data/NC_2021_Bokan_Bao_data.xlsx") %>% 
#   select(-starts_with("KI"))
data <- readxl::read_excel(snakemake@input[[1]]) %>% 
  select(-starts_with("KI"))

expr_mat <- data %>% 
  select(-IUPAC) %>% 
  column_to_rownames("glycan") %>% 
  t() %>% 
  scale() %>% 
  t()

network <- data %>% 
  select(glycan, IUPAC) %>% 
  mutate(
    FUT8 = map_int(IUPAC, n_core_fuc),
    MGAT1 = map_int(IUPAC, ~ n_glycan_type(.x) != "highmannose"),
    MGAT2 = map_int(IUPAC, ~ n_antennae(.x) > 1),
    MGAT2 = if_else(is.na(MGAT2), 0, MGAT2),
    MGAT4 = as.integer(have_motif(IUPAC, "GlcNAc(b1-4)Man(a1-")),
    MGAT5 = as.integer(have_motif(IUPAC, "GlcNAc(b1-6)Man(a1-")),
    ST3GAL = count_motif(IUPAC, "Neu5Ac(a2-"),
    B4GALT = count_motif(IUPAC, "Gal(b1-") / count_motif(IUPAC, "GlcNAc(b1-?)Man(a1-"),
    B4GALT = if_else(is.na(B4GALT), -1, B4GALT),
    B3GNT = count_motif(IUPAC, "GlcNAc(b1-3)Gal(b1-"),
  ) %>% 
  select(-IUPAC) %>% 
  pivot_longer(-glycan, names_to = "source", values_to = "mor") %>% 
  rename(target = glycan)

res <- run_ulm(expr_mat, network) %>% 
  group_by(source) %>% 
  mutate(score = as.double(scale(score))) %>% 
  ungroup() %>% 
  select(gt = source, sample = condition, score) %>% 
  pivot_wider(names_from = sample, values_from = score)
res <- res %>% 
  mutate(across(where(is.numeric), ~ .x - res$WT)) %>% 
  pivot_longer(-gt, names_to = "sample", values_to = "score")

heatmap <- ggplot(res, aes(sample, gt, fill = score)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# tgutil::ggpreview(heatmap, width = 6, height = 4)

ggsave(snakemake@output[[1]], heatmap, width = 6, height = 4)