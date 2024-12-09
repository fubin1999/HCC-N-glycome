library(tidyverse)


# data <- read_tsv("data/glycoproteome.list") %>%
#   janitor::clean_names()

data <- read_tsv(snakemake@input[[1]]) %>%
  janitor::clean_names()

groups <- tribble(
  ~raw_name, ~ratio_channel, ~group, ~sample,
  "20241017-LXJ-MIX1", "126_130C", "HC", "GPS1",
  "20241017-LXJ-MIX1", "127N_130C", "CHB", "GPS2",
  "20241017-LXJ-MIX1", "127C_130C", "LC", "GPS3",
  "20241017-LXJ-MIX1", "128N_130C", "HCC", "GPS4",
  "20241017-LXJ-MIX1", "128C_130C", "HC", "GPS5",
  "20241017-LXJ-MIX1", "129N_130C", "CHB", "GPS6",
  "20241017-LXJ-MIX1", "129C_130C", "LC", "GPS7",
  "20241017-LXJ-MIX1", "130N_130C", "HCC", "GPS8",
  "20241017-LXJ-MIX2", "126_130C", "HC", "GPS9",
  "20241017-LXJ-MIX2", "127N_130C", "CHB", "GPS10",
  "20241017-LXJ-MIX2", "127C_130C", "LC", "GPS11",
  "20241017-LXJ-MIX2", "128N_130C", "HCC", "GPS12",
  "20241017-LXJ-MIX2", "128C_130C", "HC", "GPS13",
  "20241017-LXJ-MIX2", "129N_130C", "CHB", "GPS14",
  "20241017-LXJ-MIX2", "129C_130C", "LC", "GPS15",
  "20241017-LXJ-MIX2", "130N_130C", "HCC", "GPS16",
)

commonly_identified <- data %>%
  distinct(raw_name, peptide, gly_site, glycan_composition) %>%
  summarise(n = n(), .by = c(peptide, gly_site, glycan_composition)) %>%
  filter(n == 2) %>%
  select(-n)

prepared <- data %>%
  # 0. Keep only the glycopeptides that are commonly identified in both replicates.
  semi_join(commonly_identified, by = c("peptide", "gly_site", "glycan_composition")) %>%
  # 1. Calculate log2 ratio of each report ion to the 130C report ion.
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
  ) %>%
  filter(report_ion_130C != 0) %>%
  mutate(
    ratio_126_130C = report_ion_126 / report_ion_130C,
    ratio_127N_130C = report_ion_127N / report_ion_130C,
    ratio_127C_130C = report_ion_127C / report_ion_130C,
    ratio_128N_130C = report_ion_128N / report_ion_130C,
    ratio_128C_130C = report_ion_128C / report_ion_130C,
    ratio_129N_130C = report_ion_129N / report_ion_130C,
    ratio_129C_130C = report_ion_129C / report_ion_130C,
    ratio_130N_130C = report_ion_130N / report_ion_130C,
  ) %>%
  select(
    raw_name, peptide, proteins, genes, gly_site, glycan_composition,
    ratio_126_130C, ratio_127N_130C, ratio_127C_130C, ratio_128N_130C,
    ratio_128C_130C, ratio_129N_130C, ratio_129C_130C, ratio_130N_130C
  ) %>%
  pivot_longer(
    cols = starts_with("ratio"),
    names_to = "ratio_channel",
    values_to = "ratio",
    names_prefix = "ratio_"
  ) %>%
  mutate(ratio = if_else(ratio == 0, NA_real_, ratio)) %>%
  mutate(log2_ratio = log2(ratio), .keep = "unused") %>%
  # 2. Take the median of each glycopeptide.
  summarise(
    log2_ratio = median(log2_ratio, na.rm = TRUE),
    .by = c(raw_name, ratio_channel, peptide, gly_site, glycan_composition, proteins, genes)
  ) %>%
  # 3. Median normalization.
  mutate(log2_ratio = log2_ratio - median(log2_ratio, na.rm = TRUE), .by = raw_name) %>%
  # 4. Add sample and group information.
  left_join(groups, by = c("raw_name", "ratio_channel")) %>%
  select(-raw_name, -ratio_channel) %>%
  relocate(sample, group, everything()) %>%
  # 5. Reformat some columns.
  mutate(
    glycan_composition = str_remove_all(glycan_composition, "[\\(\\)]"),
    glycan_composition = str_replace_all(glycan_composition, "A", "S"),
    proteins = str_remove(proteins, ";$"),
    genes = str_remove(genes, ";$"),
  )

write_csv(prepared, snakemake@output[[1]])