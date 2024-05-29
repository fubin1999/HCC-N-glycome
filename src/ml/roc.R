source("renv/activate.R")

library(tidyverse)
library(pROC)
library(patchwork)

# Read data-----
# predictions <- read_csv("results/data/ml/predictions.csv")
# groups <- read_csv("results/data/prepared/groups.csv")
# clinical <- read_csv("results/data/prepared/clinical.csv") |>
#   select(sample, AFP)

predictions <- read_csv(snakemake@input[["predictions"]])
groups <- read_csv(snakemake@input[["groups"]])
clinical <- read_csv(snakemake@input[["clinical"]]) |>
  select(sample, AFP)

predictions <- predictions |> 
  left_join(groups, by = "sample") |> 
  left_join(clinical, by = "sample")

# All data-----
draw_roc_curve <- function (roc_list) {
  ggroc(roc_list) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey") +
    labs(color = "") +
    theme_bw() +
    theme(
      legend.position = c(0.65, 0.25),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.background = element_blank(),
      title = element_text(size = 10)
    ) +
    scale_color_manual(values = c("#CC5F5A", "#7A848D", "#A2AFA6"))
}

global_list <- list(
  `HCC Fusion` = roc(predictions, response = "target", predictor = "probability", ci = TRUE),
  `AFP` = roc(predictions, response = "target", predictor = "AFP", ci = TRUE)
)
global_p <- draw_roc_curve(global_list) +
  ggtitle("Test Set: Control/HCC")

HC_list <- list(
  `HCC Fusion` = roc(predictions |> filter(group %in% c("HC", "HCC")), response = "target", predictor = "probability", ci = TRUE),
  `AFP` = roc(predictions |> filter(group %in% c("HC", "HCC")), response = "target", predictor = "AFP", ci = TRUE)
)
HC_p <- draw_roc_curve(HC_list) +
  ggtitle("Test Set: HC/HCC")

MC_list <- list(
  `HCC Fusion` = roc(predictions |> filter(group %in% c("CHB", "HCC")), response = "target", predictor = "probability", ci = TRUE),
  `AFP` = roc(predictions |> filter(group %in% c("CHB", "HCC")), response = "target", predictor = "AFP", ci = TRUE)
)
MC_p <- draw_roc_curve(MC_list) +
  ggtitle("Test Set: CHB/HCC")

YC_list <- list(
  `HCC Fusion` = roc(predictions |> filter(group %in% c("LC", "HCC")), response = "target", predictor = "probability", ci = TRUE),
  `AFP` = roc(predictions |> filter(group %in% c("LC", "HCC")), response = "target", predictor = "AFP", ci = TRUE)
)
YC_p <- draw_roc_curve(YC_list) +
  ggtitle("Test Set: LC/HCC")

# AFP negative samples-----
AFP_neg_samples <- predictions |> filter(AFP < 20)

AC_global_list <- list(
  `HCC Fusion` = roc(AFP_neg_samples, response = "target", predictor = "probability", ci = TRUE)
)
AN_global_p <- draw_roc_curve(AC_global_list) +
  ggtitle("Test Set: Control/HCC, AFP(-)")

AN_HC_list <- list(
  `HCC Fusion` = roc(AFP_neg_samples |> filter(group %in% c("HC", "HCC")), response = "target", predictor = "probability", ci = TRUE)
)
AN_HC_p <- draw_roc_curve(AN_HC_list) +
  ggtitle("Test Set: HC/HCC, AFP(-)")

AN_MC_list <- list(
  `HCC Fusion` = roc(AFP_neg_samples |> filter(group %in% c("CHB", "HCC")), response = "target", predictor = "probability", ci = TRUE)
)
AN_MC_p <- draw_roc_curve(AN_MC_list) +
  ggtitle("Test Set: CHB/HCC, AFP(-)")

AN_YC_list <- list(
  `HCC Fusion` = roc(AFP_neg_samples |> filter(group %in% c("LC", "HCC")), response = "target", predictor = "probability", ci = TRUE)
)
AN_YC_p <- draw_roc_curve(AN_YC_list) +
  ggtitle("Test Set: LC/HCC, AFP(-)")

final_p <- (global_p | HC_p | MC_p | YC_p) /
  (AN_global_p | AN_HC_p | AN_MC_p | AN_YC_p)
# tgutil::ggpreview(width = 12, height = 6)
ggsave(snakemake@output[[1]], final_p, width = 12, height = 6)

# Save results-----
auc_df <- tribble(
  ~data, ~between, ~predictor, ~AUC, ~ci_lower, ~ci_upper,
  "all", "Control/HCC", "HCC Fusion", global_list[["HCC Fusion"]]$auc[1], global_list[["HCC Fusion"]]$ci[1], global_list[["HCC Fusion"]]$ci[3],
  "all", "Control/HCC", "AFP", global_list[["AFP"]]$auc[1], global_list[["AFP"]]$ci[1], global_list[["AFP"]]$ci[3],
  "all", "HC/HCC", "HCC Fusion", HC_list[["HCC Fusion"]]$auc[1], HC_list[["HCC Fusion"]]$ci[1], HC_list[["HCC Fusion"]]$ci[3],
  "all", "HC/HCC", "AFP", HC_list[["AFP"]]$auc[1], HC_list[["AFP"]]$ci[1], HC_list[["AFP"]]$ci[3],
  "all", "CHB/HCC", "HCC Fusion", MC_list[["HCC Fusion"]]$auc[1], MC_list[["HCC Fusion"]]$ci[1], MC_list[["HCC Fusion"]]$ci[3],
  "all", "CHB/HCC", "AFP", MC_list[["AFP"]]$auc[1], MC_list[["AFP"]]$ci[1], MC_list[["AFP"]]$ci[3],
  "all", "LC/HCC", "HCC Fusion", YC_list[["HCC Fusion"]]$auc[1], YC_list[["HCC Fusion"]]$ci[1], YC_list[["HCC Fusion"]]$ci[3],
  "all", "LC/HCC", "AFP", YC_list[["AFP"]]$auc[1], YC_list[["AFP"]]$ci[1], YC_list[["AFP"]]$ci[3],
  "AFP-", "Control/HCC", "HCC Fusion", AC_global_list[["HCC Fusion"]]$auc[1], AC_global_list[["HCC Fusion"]]$ci[1], AC_global_list[["HCC Fusion"]]$ci[3],
  "AFP-", "HC/HCC", "HCC Fusion", AN_HC_list[["HCC Fusion"]]$auc[1], AN_HC_list[["HCC Fusion"]]$ci[1], AN_HC_list[["HCC Fusion"]]$ci[3],
  "AFP-", "CHB/HCC", "HCC Fusion", AN_MC_list[["HCC Fusion"]]$auc[1], AN_MC_list[["HCC Fusion"]]$ci[1], AN_MC_list[["HCC Fusion"]]$ci[3],
  "AFP-", "LC/HCC", "HCC Fusion", AN_YC_list[["HCC Fusion"]]$auc[1], AN_YC_list[["HCC Fusion"]]$ci[1], AN_YC_list[["HCC Fusion"]]$ci[3]
)
write_csv(auc_df, snakemake@output[[2]])
