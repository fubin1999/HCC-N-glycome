library(tidyverse)

# raw_clinical <- read_csv("data/clinical.csv")
# plates <- read_csv("data/plates.csv")
# groups <- read_csv("results/data/prepared/unfiltered_groups.csv")

raw_clinical <- read_csv(snakemake@input[["clinical"]])
plates <- read_csv(snakemake@input[["plates"]])
groups <- read_csv(snakemake@input[["groups"]])

# Convert sample names and keep only samples that have abundance data
clinical <- raw_clinical |> 
  left_join(plates |> select(raw_sample, sample), by = "raw_sample") |>
  filter(!is.na(sample)) |>
  select(sample, everything()) |>
  select(-raw_sample) |> 
  mutate(sample_no = parse_number(sample)) |>
  arrange(sample_no) |> 
  select(-sample_no)

# Clean up clinical data
clean_clinical <- clinical %>%
  distinct(sample, .keep_all = TRUE) %>%
  mutate(across(-c(sample, sex, child_pugh, TNM), ~ parse_number(as.character(.x)))) %>%
  # Convert some clinical variables from numeric to postive/negative factors
  mutate(
    HBSAG = case_when(
      is.na(HBSAG) ~ "-",
      HBSAG < 0.5 ~ "-",
      .default = "+"
    ),
    HBEAG = case_when(
      is.na(HBEAG) ~ "-",
      HBEAG < 0.5 ~ "-",
      .default = "+"
    ),
    HBEAB = case_when(
      is.na(HBEAB) ~ "-",
      HBEAB < 0.2 ~ "-",
      .default = "+"
    ),
    HBCAB = case_when(
      is.na(HBCAB) ~ "-",
      HBCAB < 0.9 ~ "-",
      .default = "+"
    ),
    HCV = case_when(
      is.na(HCV) ~ "-",
      HCV < 1000 ~ "-",
      .default = "+"
    ),
  ) %>%
  # Extract TNM stage
  mutate(
    T_stage = str_extract(TNM, "T([0-4X])", group = 1),
    T_stage = if_else(T_stage == "X", NA_character_, T_stage),
    N_stage = str_extract(TNM, "N([0-3X])", group = 1),
    N_stage = if_else(N_stage == "X", NA_character_, N_stage),
    M_stage = str_extract(TNM, "M([0-1X])", group = 1),
    M_stage = if_else(M_stage == "X", NA_character_, M_stage)
  ) %>%
  mutate(
    TNM_stage = case_when(
      (N_stage %in% c("1", "2", "3")) | (M_stage == "1") ~ "IV",
      T_stage %in% c("3", "4") ~ "III",
      T_stage == "2" ~ "II",
      T_stage == "1" ~ "I",
      .default = NA_character_
    )
  ) %>%
  # Here we remove some variables that are not useful for our analysis.
  # TNM, T_stage, N_stage, M_stage are used to calculate TNM_stage, so no longer needed.
  # CEA, CA199, and HCV are not related to the research question.
  select(-c(TNM, T_stage, N_stage, M_stage, CEA, CA199, HCV))

# Impute missing values using KNN
imputed_clinical <- clean_clinical %>%
  mutate(AFP = replace_na(AFP, 0)) %>%
  VIM::kNN(imp_var = FALSE, variable = c("sex", "age", "AST", "ALT", "GGT", "ALB", "TBIL", "TP", "child_pugh")) %>%
  as_tibble()

# Impute TNM stage for HCC samples
HCC_imputed_clinical <- imputed_clinical %>%
  left_join(groups, by = "sample") %>%
  filter(group == "HCC") %>%
  select(-group)

non_HCC_imputed_clinical <- imputed_clinical %>%
  left_join(groups, by = "sample") %>%
  filter(group != "HCC") %>%
  select(-group)

HCC_imputed_clinical_final <- HCC_imputed_clinical %>%
  VIM::kNN(imp_var = FALSE, variable = "TNM_stage") %>%
  as_tibble()

imputed_clinical_final <- bind_rows(HCC_imputed_clinical_final, non_HCC_imputed_clinical)

# Calculate AAR and ALBI
final_clinical <- imputed_clinical_final |>
  mutate(AAR = AST / ALT) %>%
  mutate(
    ALBI = log10(TBIL) * 0.66 - ALB * 0.085,
    ALBI_stage = case_when(
      ALBI <= -2.60 ~ "I",
      ALBI > -2.60 & ALBI <= -1.39 ~ "II",
      ALBI > -1.39 ~ "III",
    )
  ) %>%
  select(-ALBI)

write_csv(final_clinical, snakemake@output[[1]])
