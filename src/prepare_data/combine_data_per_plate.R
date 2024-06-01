library(tidyverse)

convert <- function(comp) {
  comp = str_replace(comp, "\\[.*?\\]", "")
  comp = str_replace(comp, "Hex\\((\\d+)\\)", "H\\1")
  comp = str_replace(comp, "HexNAc\\((\\d+)\\)", "N\\1")
  comp = str_replace(comp, "dHex\\((\\d+)\\)", "F\\1")
  comp = str_replace(comp, "NeuAc\\((\\d+)\\)", "S\\1")
  return(comp)
}

files <- unlist(snakemake@input)
combined <- read_csv(files) |> 
  complete(sample, glycan, fill = list(value = NA)) |>
  mutate(glycan = convert(glycan))
wide_data <- combined |>
  pivot_wider(names_from = glycan, values_from = value)
write_csv(wide_data, snakemake@output[[1]])