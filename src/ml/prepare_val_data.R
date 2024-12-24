library(tidyverse)
library(sva)

# Load the data-----
# raw_data <- read_csv("results/data/glyhunter_results_GD2/summary_area.csv")
# sample_info <- read_csv("data/cohort_GD2/maldi_pos_and_group.csv")
# train_data <- read_csv("results/data/ml/train_data.csv")

raw_data <- read_csv(snakemake@input[[1]])
sample_info <- read_csv(snakemake@input[[2]])
train_data <- read_csv(snakemake@input[[3]])

sample_info <- sample_info %>% 
  select(sample, maldi_pos, group)

# Prepare the data-----
convert <- function(comp) {
  comp <- str_replace(comp, "\\[.*?\\]", "")
  comp <- str_replace(comp, "Hex\\((\\d+)\\)", "H\\1")
  comp <- str_replace(comp, "HexNAc\\((\\d+)\\)", "N\\1")
  comp <- str_replace(comp, "dHex\\((\\d+)\\)", "F\\1")
  comp <- str_replace(comp, "NeuAc\\((\\d+)\\)", "S\\1")
  return(comp)
}

selected <- c('H3N4F1', 'H3N5F1', 'H4N3S1', 'H4N4', 'H4N4S1', 'H4N4F1', 'H4N5F1', 'H5N2', 'H5N4',
             'H5N4S1', 'H5N4S2', 'H5N4F1', 'H5N4F1S1', 'H5N4F1S2', 'H5N5S1', 'H5N5F1', 'H5N5F1S1',
             'H5N5F1S2', 'H6N2', 'H6N5S2', 'H6N5S3', 'H6N5F1S2', 'H6N5F1S3', 'H7N2', 'H8N2', 'H9N2')

prepared <- raw_data %>%
  mutate(glycan = convert(glycan)) %>%
  pivot_longer(-glycan, names_to = "sample", values_to = "value") %>%
  filter(glycan %in% selected) %>%
  mutate(value = value / sum(value, na.rm = TRUE) * 100, .by = sample) %>%
  separate_wider_delim(sample, "_", names = c(NA, NA, "maldi_pos", NA)) %>%
  left_join(sample_info, by = "maldi_pos") %>%
  select(-maldi_pos) %>%
  pivot_wider(names_from = glycan, values_from = value) %>%
  filter(group != "QC") %>%
  mutate(group = case_match(group, "H" ~ "HC", "M" ~ "CHB", "Y" ~ "LC", "C" ~ "HCC"))

batch_corrected <- local({
  all_data <- bind_rows(
    train_data %>%
      relocate(sample, group, all_of(selected)) %>%
      mutate(source = "train", .before = 1),
    prepared %>%
      relocate(sample, group, all_of(selected)) %>%
      mutate(source = "test", .before = 1)
  )
  expr_mat <- all_data %>%
    select(-source, -group) %>%
    column_to_rownames("sample") %>%
    as.matrix() %>%
    t() %>%
    log()
  sample_info <- all_data %>%
    distinct(sample, group, source) %>%
    column_to_rownames("sample") %>%
    .[colnames(expr_mat), c("source", "group")] %>%
    mutate(source = factor(source), group = factor(group))
  mod_combat <- model.matrix(~group, data = sample_info)
  combat_expr_mat <- ComBat(dat = expr_mat, batch = sample_info$source, mod = mod_combat, ref.batch = "train")
  combat_all_data <- combat_expr_mat %>%
    exp() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    as_tibble() %>%
    left_join(all_data %>% distinct(sample, group, source), by = "sample")
  combat_all_data %>%
    filter(source == "test") %>%
    select(-source) %>%
    relocate(sample, group, all_of(selected))
})

# Save the data-----
write_csv(batch_corrected, snakemake@output[[1]])