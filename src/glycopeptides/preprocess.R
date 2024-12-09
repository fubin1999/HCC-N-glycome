library(tidyverse)


data <- read_tsv("data/glycoproteome.list") %>%
  janitor::clean_names()

data <- read_tsv(snakemake@input[[1]]) %>%
  janitor::clean_names()

groups <- tribble(
  ~raw_name, ~ratio, ~group, ~sample,
  "20241017-LXJ-MIX1", "ratio_127N_126", "HC", "GPS1",
  "20241017-LXJ-MIX1", "ratio_127C_126", "CHB", "GPS2",
  "20241017-LXJ-MIX1", "ratio_128N_126", "LC", "GPS3",
  "20241017-LXJ-MIX1", "ratio_128C_126", "HCC", "GPS4",
  "20241017-LXJ-MIX1", "ratio_129N_126", "HC", "GPS5",
  "20241017-LXJ-MIX1", "ratio_129C_126", "CHB", "GPS6",
  "20241017-LXJ-MIX1", "ratio_130N_126", "LC", "GPS7",
  "20241017-LXJ-MIX1", "ratio_130C_126", "HCC", "GPS8",
  "20241017-LXJ-MIX2", "ratio_127N_126", "HC", "GPS9",
  "20241017-LXJ-MIX2", "ratio_127C_126", "CHB", "GPS10",
  "20241017-LXJ-MIX2", "ratio_128N_126", "LC", "GPS11",
  "20241017-LXJ-MIX2", "ratio_128C_126", "HCC", "GPS12",
  "20241017-LXJ-MIX2", "ratio_129N_126", "HC", "GPS13",
  "20241017-LXJ-MIX2", "ratio_129C_126", "CHB", "GPS14",
  "20241017-LXJ-MIX2", "ratio_130N_126", "LC", "GPS15",
  "20241017-LXJ-MIX2", "ratio_130C_126", "HCC", "GPS16",
)

prepared <- data %>%
  rename(
    report_ion_126 = report_ion_126_1277,
    report_ion_127N = report_ion_127_1248,
    report_ion_127C = report_ion_127_1311,
    report_ion_128N = report_ion_128_1281,
    report_ion_128C = report_ion_128_1344,
    report_ion_129N = report_ion_129_1315,
    report_ion_129C = report_ion_129_1378,
    report_ion_130N = report_ion_130_1348,
    report_ion_130C = report_ion_130_1411,
    report_ion_131N = report_ion_131_1382,
    report_ion_131C = report_ion_131_1445
  ) %>%
  filter(report_ion_126 != 0) %>%
  mutate(
    ratio_127N_126 = report_ion_127N / report_ion_126,
    ratio_127C_126 = report_ion_127C / report_ion_126,
    ratio_128N_126 = report_ion_128N / report_ion_126,
    ratio_128C_126 = report_ion_128C / report_ion_126,
    ratio_129N_126 = report_ion_129N / report_ion_126,
    ratio_129C_126 = report_ion_129C / report_ion_126,
    ratio_130N_126 = report_ion_130N / report_ion_126,
    ratio_130C_126 = report_ion_130C / report_ion_126,
  ) %>%
  mutate(across(starts_with("ratio"), ~ log2(.))) %>%
  mutate(across(starts_with("ratio"), ~ . - median(.))) %>%
  summarise(
    across(starts_with("ratio"), ~ median(.)),
    .by = c(raw_name, peptide, proteins, genes, gly_site, glycan_composition)
  ) %>%
  pivot_longer(starts_with("ratio"), names_to = "ratio", values_to = "log2_ratio") %>%
  left_join(groups, by = c("raw_name", "ratio")) %>%
  select(-raw_name, -ratio) %>%
  mutate(
    proteins = str_remove(proteins, ";$"),
    genes = str_remove(genes, ";$"),
    glycan_composition = str_remove_all(glycan_composition, "[\\(\\)]"),
    glycan_composition = str_replace(glycan_composition, "A", "S")
  ) %>%
  select(-group) %>%
  pivot_wider(names_from = sample, values_from = log2_ratio) %>%
  pivot_longer(starts_with("GPS"), names_to = "sample", values_to = "log2_ratio") %>%
  left_join(groups %>% select(sample, group), by = "sample")

write_csv(prepared, snakemake@output[[1]])