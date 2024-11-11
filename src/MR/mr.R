library(TwoSampleMR)
library(tidyverse)

glycan_exp_df <- read.csv("data/GWAS_cwaa053.csv")
glycan_exp_dat <- format_data(glycan_exp_df, type = "exposure")
glycan_exp_dat <- clump_data(glycan_exp_dat)
hcc_out_dat <- extract_outcome_data(snps = glycan_exp_dat$SNP, outcomes = c("ebi-a-GCST90092003", "bbj-a-158"), maf_threshold = 0.1)
dat <- harmonise_data(exposure_dat = glycan_exp_dat, outcome_dat = hcc_out_dat)
res <- mr(dat)

write_csv(res, "results/data/MR/MR_result.csv")
