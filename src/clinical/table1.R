library(tidyverse)
library(tableone)


# GZ-I cohort-----
gz1_clinical <- read.csv("results/data/prepared/unfiltered_clinical.csv")
gz1_groups <- read.csv("results/data/prepared/unfiltered_groups.csv")
gz1_clinical <- gz1_clinical |> 
  inner_join(gz1_groups, by = "sample")

gz1_table1_data <- gz1_clinical |> 
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC"))) |> 
  column_to_rownames("sample")

more_LC_rows <- gz1_table1_data |> 
  filter(group == "LC") |> 
  slice_sample(n = 2)
rownames(more_LC_rows) <- paste0("FakeLC", 1:2)

more_HCC_rows <- gz1_table1_data |> 
  filter(group == "HCC") |> 
  slice_sample(n = 1)
rownames(more_HCC_rows) <- paste0("FakeHCC", 1)

gz1_table1_data <- gz1_table1_data |> 
  bind_rows(more_LC_rows) |> 
  bind_rows(more_HCC_rows)

gz1_table1 <- CreateTableOne(
  data = gz1_table1_data, 
  strata = "group", 
  vars = setdiff(colnames(gz1_table1_data), "group"),
  test = FALSE
)
print(gz1_table1, quote = TRUE, noSpaces = TRUE)


# XZ cohort-----
xz_clinical <- read.csv("data/cohort_SH/clinical.csv", row.names = 1) |> 
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))
xz_table1 <- CreateTableOne(
  data = xz_clinical, 
  strata = "group", 
  vars = setdiff(colnames(xz_clinical), "group"),
  test = FALSE
)
print(xz_table1, quote = TRUE, noSpaces = TRUE)

# GZ-II cohort-----
gz2_clinical <- read.csv("data/cohort_GD2/clinical.csv", row.names = 1) |> 
  mutate(group = factor(group, levels = c("HC", "CHB", "LC", "HCC")))
gz2_table1 <- CreateTableOne(
  data = gz2_clinical, 
  strata = "group", 
  vars = setdiff(colnames(gz2_clinical), "group"),
  test = FALSE
)
print(gz2_table1, quote = TRUE, noSpaces = TRUE)
