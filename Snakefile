PREPARED_DIR = "results/data/prepared/"
RAW_ABUNDANCE = PREPARED_DIR + "raw_abundance.csv"
PROCESSED_ABUNDANCE = PREPARED_DIR + "processed_abundance.csv"
ALL_DERIVED_TRAITS = PREPARED_DIR + "derived_traits.csv"
FILTERED_DERIVED_TRAITS = PREPARED_DIR + "filtered_derived_traits.csv"
META_PROPERTIES = PREPARED_DIR + "meta_properties.csv"
GROUPS = PREPARED_DIR + "groups.csv"
CLINICAL = PREPARED_DIR + "clinical.csv"

rule all:
    input:
        # ===== Others =====
        "results/data/prepared/raw_abundance_full.csv",
        "results/data/clinical/AFP_cutoff.csv",
        "results/figures/clinical/AFP_cutoff.pdf",
        "results/figures/clinical/clinical_heatmap.pdf",
        "results/figures/SNFG/complete_svg/",
        "results/figures/SNFG/compact_svg/",
        "results/figures/SNFG/complete_pdf/",
        "results/figures/SNFG/compact_pdf/",

        # ===== Data Quality Figures =====
        "results/figures/data_quality/batch_effect_pca.pdf",
        "results/figures/data_quality/glycan_count_venn_all.pdf",
        "results/figures/data_quality/glycan_count_venn_with_structures.pdf",
        "results/figures/data_quality/glycan_count_venn_confident.pdf",
        "results/figures/data_quality/glycan_property_heatmap.pdf",
        "results/figures/data_quality/glycan_property_barplots.pdf",
        "results/figures/data_quality/glycan_property_boxplots.pdf",
        "results/figures/data_quality/QC_corrplot.pdf",

        # ===== Differential Analysis Data =====
        "results/data/diff_analysis/detect_rate_diff.csv",
        "results/data/diff_analysis/glycan_ancova.csv",
        "results/data/diff_analysis/glycan_post_hoc.csv",
        "results/data/diff_analysis/glycan_fold_change.csv",
        "results/data/diff_analysis/trait_ancova.csv",
        "results/data/diff_analysis/trait_post_hoc.csv",
        "results/data/diff_analysis/trait_fold_change.csv",

        # ===== Differential Analysis Figures =====
        "results/figures/diff_analysis/detect_rate_diff.pdf",
        "results/figures/diff_analysis/glycan_diff_rose_plot.pdf",
        "results/figures/diff_analysis/glycan_violin_plots.pdf",
        "results/figures/diff_analysis/glycan_volcanos_lf_adjusted.pdf",
        "results/figures/diff_analysis/glycan_diff_bubble.pdf",
        "results/figures/diff_analysis/glycan_heatmap.pdf",
        "results/figures/diff_analysis/glycan_volcanos.pdf",
        "results/figures/diff_analysis/glycan_compare_FC.pdf",
        "results/figures/diff_analysis/glycan_compare_p.pdf",
        "results/figures/diff_analysis/glycan_confounders.pdf",
        "results/figures/diff_analysis/glycan_pca.pdf",
        "results/figures/diff_analysis/trait_confounders.pdf",
        "results/figures/diff_analysis/trait_boxplots.pdf",
        "results/figures/diff_analysis/trait_diff_bubble.pdf",
        "results/figures/diff_analysis/trait_heatmap.pdf",

        # ===== Glycan Coexpression Module Data =====
        "results/figures/glycan_coexpr/cc_result/",
        "results/data/glycan_coexpr/glycan_clusters.csv",
        "results/data/glycan_coexpr/eigen_glycans.csv",
        "results/data/glycan_coexpr/cluster_ancova.csv",
        "results/data/glycan_coexpr/cluster_post_hoc.csv",

        # ===== Glycan Coexpression Module Figures =====
        "results/figures/glycan_coexpr/cluster_glycan_heatmap.pdf",
        "results/figures/glycan_coexpr/glycan_cluster_trends.pdf",
        "results/figures/glycan_coexpr/glycan_property_heatmap.pdf",
        "results/figures/glycan_coexpr/cluster_corrplot.pdf",
        "results/figures/glycan_coexpr/cor_inter_intra_GCM.pdf",

        # ===== Correlation with Liver Function Data =====
        "results/data/cor_with_liver_function/global_cor_result_with_liver_functions.csv",
        "results/data/cor_with_liver_function/grouped_cor_result_with_liver_functions.csv",
        "results/data/cor_with_liver_function/global_ttest_result_with_liver_functions.csv",
        "results/data/cor_with_liver_function/grouped_ttest_result_with_liver_functions.csv",
        "results/data/cor_with_liver_function/ALBI_model_scores.csv",
        "results/data/cor_with_liver_function/ALBI_model_preds.csv",
        "results/data/cor_with_liver_function/ALBI_score_model_metrics.csv",

        # ===== Correlation with Liver Function Figures =====
        "results/figures/cor_with_liver_function/trait_clinical_subtype_boxplots.pdf",
        "results/figures/cor_with_liver_function/glycan_liver_function_heatmap.pdf",
        "results/figures/cor_with_liver_function/cor_coef_per_group.pdf",
        "results/figures/cor_with_liver_function/ALBI_model_ROC.pdf",
        "results/figures/cor_with_liver_function/ALBI_model_confusion_matrix.pdf",
        "results/figures/cor_with_liver_function/ALBI_model_confusion_matrix_grouped.pdf",
        "results/figures/cor_with_liver_function/ALBI_model_prob_ridge_plot.pdf",
        "results/figures/cor_with_liver_function/ALBI_score_model_prediction.pdf",

        # ===== Molecular Subtypes Data =====
        "results/data/subtypes/consensus_cluster_result.csv",
        "results/data/subtypes/subtype_glycan_anova.csv",
        "results/data/subtypes/subtype_glycan_post_hoc.csv",
        "results/data/subtypes/subtype_with_continous_clinical_kruskal_result.csv",
        "results/data/subtypes/subtype_with_continous_clinical_post_hoc_result.csv",
        "results/data/subtypes/subtype_with_categoric_clinical_fisher_result.csv",
        "results/data/subtypes/subtype_with_categoric_clinical_post_hoc_result.csv",
        "results/data/subtypes/subtype_coexp_module_anova.csv",
        "results/data/subtypes/subtype_coexp_module_post_hoc.csv",
        "results/data/subtypes/subtype_with_other_groups_glycan_ttest.csv",
        "results/data/subtypes/subtype_with_other_groups_clinical_wilcox.csv",

        # ===== Molecular Subtypes Figures =====
        "results/figures/subtypes/cc_result/",
        "results/figures/subtypes/subtype_heatmap.pdf",
        "results/figures/subtypes/subtype_pca.pdf",
        "results/figures/subtypes/subtype_coexp_module_heatmap.pdf",
        "results/figures/subtypes/subtype_continuous_clinical_boxplots.pdf",
        "results/figures/subtypes/hypergeometric_test_heatmap.pdf",
        "results/figures/subtypes/fisher_test_barplot.pdf",
        "results/figures/subtypes/compare_with_other_group_heatmap.pdf",
        "results/figures/subtypes/compare_with_other_group_glycan_diff_heatmap.pdf",
        "results/figures/subtypes/compare_with_other_group_clinical_diff_heatmap.pdf",

        # ===== TCGA Data =====
        "results/data/TCGA/dea_results.csv",
        "results/data/TCGA/consensus_cluster_result.csv",
        "results/data/TCGA/cluster_dea_results.csv",
        "results/data/TCGA/cox.csv",

        # ===== TCGA Figures =====
        "results/figures/TCGA/volcano.pdf",
        "results/figures/TCGA/MA_plot.pdf",
        "results/figures/TCGA/heatmap.pdf",
        "results/figures/TCGA/consensus_cluster/",
        "results/figures/TCGA/cluster_dea_heatmap.pdf",
        "results/figures/TCGA/clinical_heatmap.pdf",
        "results/figures/TCGA/cluster_KM.pdf",
        "results/figures/TCGA/single_gene_KM/",

# ==================== Prepare Data ====================
rule build_db:
    # Convert the serum glycan CSV file into a byonic database for GlyHunter.
    input:
        "data/human_serum_glycans.csv"
    output:
        "results/data/glycan_db.byonic"
    script:
        "src/prepare_data/make_db.py"

rule run_glyhunter:
    # Run GlyHunter on each mass list file.
    input:
        mass_list="data/mass_list/plate{no}.xlsx",
        db="results/data/glycan_db.byonic"
    output:
        directory("results/data/glyhunter_results/plate{no}/")
    shell:
        "glyhunter run {input.mass_list} -d {input.db} -o {output}"

rule run_glyhunter_default_db:
    # Run GlyHunter with default serum N-glycan database.
    input:
        "data/mass_list/plate{no}.xlsx"
    output:
        directory("results/data/glyhunter_results_full/plate{no}/")
    shell:
        "glyhunter run {input} -o {output}"

rule assign_maldi_pos:
    # Match the MALDI positions of GlyHunter results to samples.
    input:
        "data/MALDI_positions.csv",
        "results/data/glyhunter_results/plate{no}/"
    output:
        "results/data/data_per_plate/plate{no}.csv"
    script:
        "src/prepare_data/assign_maldi_pos.R"

rule assign_maldi_pos_full:
    # Match the MALDI positions of GlyHunter results to samples.
    input:
        "data/MALDI_positions.csv",
        "results/data/glyhunter_results_full/plate{no}/"
    output:
        "results/data/data_per_plate_full/plate{no}.csv"
    script:
        "src/prepare_data/assign_maldi_pos.R"

rule combine_data_per_plate:
    # Combine all GlyHunter results into a single file.
    input:
        expand("results/data/data_per_plate/plate{no}.csv", no=range(1, 9))
    output:
        RAW_ABUNDANCE
    script:
        "src/prepare_data/combine_data_per_plate.R"

rule combine_data_per_plate_full:
    # Combine all GlyHunter results into a single file.
    input:
        expand("results/data/data_per_plate_full/plate{no}.csv", no=range(1, 9))
    output:
        "results/data/prepared/raw_abundance_full.csv"
    script:
        "src/prepare_data/combine_data_per_plate.R"

rule prepare_groups:
    # Prepare the groups.
    input:
        plates="data/plates.csv"
    output:
        PREPARED_DIR + "unfiltered_groups.csv"
    script:
        "src/prepare_data/prepare_groups.R"

rule prepare_clinical:
    # Prepare the clinical information.
    input:
        clinical="data/clinical.csv",
        plates="data/plates.csv",
        groups=PREPARED_DIR + "unfiltered_groups.csv"
    output:
        PREPARED_DIR + "unfiltered_clinical.csv"
    script:
        "src/prepare_data/prepare_clinical.R"

rule preprocess:
    # Filter glycan, impute missing values, and normalize.
    input:
        RAW_ABUNDANCE,
        PREPARED_DIR + "unfiltered_groups.csv",
        PREPARED_DIR + "unfiltered_clinical.csv",
        "data/plates.csv"
    output:
        PROCESSED_ABUNDANCE,
        GROUPS,
        CLINICAL
    script:
        "src/prepare_data/preprocess.R"

rule calculate_derived_traits:
    # Calculate derived traits using GlyTrait
    input:
        PROCESSED_ABUNDANCE,
        "data/human_serum_glycans.csv",
        "src/prepare_data/struc_builtin_formulas.txt"
    output:
        filtered_traits=FILTERED_DERIVED_TRAITS,
        all_traits=ALL_DERIVED_TRAITS,
        mp_table=META_PROPERTIES
    script:
        "src/prepare_data/derive_traits.py"


# ==================== Data Quality ====================
rule batch_effect_pca:
    # Draw PCA plot to check batch effect.
    input:
        PROCESSED_ABUNDANCE,
        "data/plates.csv"
    output:
        "results/figures/data_quality/batch_effect_pca.pdf"
    script:
        "src/data_quality/batch_effect.R"

rule glycan_count_venn:
    # Draw venn diagrams for glycan count per group.
    input:
        "results/data/prepared/raw_abundance_full.csv",
        "results/data/prepared/raw_abundance.csv",
        "results/data/prepared/processed_abundance.csv",
        "results/data/prepared/groups.csv"
    output:
        "results/figures/data_quality/glycan_count_venn_all.pdf",
        "results/figures/data_quality/glycan_count_venn_with_structures.pdf",
        "results/figures/data_quality/glycan_count_venn_confident.pdf"
    script:
        "src/data_quality/glycan_count_venn.R"

rule glycan_property_heatmap:
    # Draw heatmap to show all glycans' properties.
    input:
        PROCESSED_ABUNDANCE,
        RAW_ABUNDANCE,
        META_PROPERTIES
    output:
        "results/figures/data_quality/glycan_property_heatmap.pdf"
    script:
        "src/data_quality/glycan_property_heatmap.R"

rule glycan_property_detailed:
    # Draw barplots and boxplots for more detailed properties.
    # Barplots: for number of glycans in each category.
    # Boxplots: for relative abundance of glycans in each category.
    input:
        PROCESSED_ABUNDANCE,
        META_PROPERTIES
    output:
        "results/figures/data_quality/glycan_property_barplots.pdf",
        "results/figures/data_quality/glycan_property_boxplots.pdf"
    script:
        "src/data_quality/glycan_property_detailed.R"

rule QC_corrplot:
    # Draw correlation plot for QC samples.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS
    output:
        "results/figures/data_quality/QC_corrplot.pdf"
    script:
        "src/data_quality/QC_corrplot.R"


# ==================== Clinical Information ====================
rule AFP_cutoff:
    # Calculate classification metrics of AFP with different thresholds.
    input:
        CLINICAL,
        GROUPS
    output:
        "results/data/clinical/AFP_cutoff.csv",
        "results/figures/clinical/AFP_cutoff.pdf"
    script:
        "src/clinical/AFP_cutoff.R"

rule clinical_heatmap:
    # Draw heatmap to show clinical information.
    input:
        CLINICAL,
        GROUPS
    output:
        "results/figures/clinical/clinical_heatmap.pdf"
    script:
        "src/clinical/clinical_heatmap.R"


# ==================== Differential Analysis ====================
rule detection_rate_diff:
    # Perform Fisher's exact test to test detection rate difference of glycans.
    input:
        RAW_ABUNDANCE,
        GROUPS
    output:
        "results/data/diff_analysis/detect_rate_diff.csv",
        "results/figures/diff_analysis/detect_rate_diff.pdf"
    script:
        "src/diff_analysis/glycan_detect_rate.R"

rule glycan_ancova:
    # Perform ANCOVA for each glycan.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        clinical=CLINICAL
    output:
        ancova="results/data/diff_analysis/glycan_ancova.csv",
        posthoc="results/data/diff_analysis/glycan_post_hoc.csv",
        ancova_lf_adjusted="results/data/diff_analysis/glycan_ancova_lf_adjusted.csv",
        posthoc_lf_adjusted="results/data/diff_analysis/glycan_post_hoc_lf_adjusted.csv"
    script:
        "src/diff_analysis/glycan_ancova.R"

rule glycan_fold_change:
    # Calculate fold changes for each glycan.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS
    output:
        "results/data/diff_analysis/glycan_fold_change.csv"
    script:
        "src/diff_analysis/glycan_fold_change.R"

rule glycan_diff_rose_plot:
    # Draw rose plot for differential glycans between each group pair.
    input:
        "results/data/diff_analysis/glycan_ancova.csv",
        "results/data/diff_analysis/glycan_post_hoc.csv"
    output:
        "results/figures/diff_analysis/glycan_diff_rose_plot.pdf"
    script:
        "src/diff_analysis/glycan_diff_rose_plot.R"

rule glycan_diff_bubble:
    # Draw bubble plot for p-values and fold changes of glycans.
    input:
        "results/data/diff_analysis/glycan_ancova.csv",
        "results/data/diff_analysis/glycan_post_hoc.csv",
        "results/data/diff_analysis/glycan_fold_change.csv",
    output:
        "results/figures/diff_analysis/glycan_diff_bubble.pdf"
    script:
        "src/diff_analysis/glycan_diff_bubble.R"

rule glycan_heatmap:
    # Draw mean heatmap for differential glycans.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS,
        "results/data/diff_analysis/glycan_ancova.csv"
    output:
        "results/figures/diff_analysis/glycan_heatmap.pdf"
    script:
        "src/diff_analysis/glycan_heatmap.R"

rule glycan_volcano:
    # Draw glycan volcano plots.
    input:
        "results/data/diff_analysis/glycan_ancova.csv",
        "results/data/diff_analysis/glycan_post_hoc.csv",
        "results/data/diff_analysis/glycan_fold_change.csv"
    output:
        "results/figures/diff_analysis/glycan_volcanos.pdf"
    script:
        "src/diff_analysis/glycan_volcano.R"

rule glycan_volcano_lf_adjusted:
    # Draw glycan volcano plots with liver function adjusted.
    input:
        "results/data/diff_analysis/glycan_ancova_lf_adjusted.csv",
        "results/data/diff_analysis/glycan_post_hoc_lf_adjusted.csv",
        "results/data/diff_analysis/glycan_fold_change.csv"
    output:
        "results/figures/diff_analysis/glycan_volcanos_lf_adjusted.pdf"
    script:
        "src/diff_analysis/glycan_volcano.R"

rule glycan_violin_plots:
    # Plot violin plots for all significant glycans.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        ancova_result="results/data/diff_analysis/glycan_ancova.csv"
    output:
        "results/figures/diff_analysis/glycan_violin_plots.pdf"
    script:
        "src/diff_analysis/glycan_violin_plots.R"

rule compare_glycan_FC:
    # Compare the fold changes of glycans between each group pair.
    input:
        "results/data/diff_analysis/glycan_fold_change.csv"
    output:
        "results/figures/diff_analysis/glycan_compare_FC.pdf"
    script:
        "src/diff_analysis/compare_FC.R"

rule compare_glycan_p:
    # Compare the p-values of glycans between each group pair.
    input:
        "results/data/diff_analysis/glycan_post_hoc.csv"
    output:
        "results/figures/diff_analysis/glycan_compare_p.pdf"
    script:
        "src/diff_analysis/compare_p.R"

rule glycan_confounders:
    # Plot dot plot for confounders' p-values.
    input:
        "results/data/diff_analysis/glycan_ancova.csv"
    output:
        "results/figures/diff_analysis/glycan_confounders.pdf"
    script:
        "src/diff_analysis/glycan_confounder_dot_plot.R"

rule glycan_pca:
    # Draw PCA plots for glycan abundance.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS
    output:
        "results/figures/diff_analysis/glycan_pca.pdf"
    script:
        "src/diff_analysis/glycan_pca.R"

rule trait_ancova:
    # Perform ANCOVA on derived traits.
    input:
        traits=FILTERED_DERIVED_TRAITS,
        groups=GROUPS,
        clinical=CLINICAL
    output:
        "results/data/diff_analysis/trait_ancova.csv",
        "results/data/diff_analysis/trait_post_hoc.csv"
    script:
        "src/diff_analysis/trait_ancova.R"

rule trait_fold_change:
    # Calculate fold changes for each derived trait.
    input:
        FILTERED_DERIVED_TRAITS,
        GROUPS
    output:
        "results/data/diff_analysis/trait_fold_change.csv"
    script:
        "src/diff_analysis/trait_fold_change.R"

rule trait_confounders:
    # Plot dot plot for confounders' p-values.
    input:
        "results/data/diff_analysis/trait_ancova.csv"
    output:
        "results/figures/diff_analysis/trait_confounders.pdf"
    script:
        "src/diff_analysis/trait_confounder_dot_plot.R"

rule trait_boxplots:
    # Draw boxplots for selected derived traits.
    input:
        ALL_DERIVED_TRAITS,
        GROUPS,
        "results/data/diff_analysis/trait_post_hoc.csv"
    output:
       "results/figures/diff_analysis/trait_boxplots.pdf"
    script:
        "src/diff_analysis/trait_boxplots.R"

rule trait_diff_bubble:
    # Draw bubble plot for differential derived traits.
    input:
        "results/data/diff_analysis/trait_post_hoc.csv"
    output:
        "results/figures/diff_analysis/trait_diff_bubble.pdf"
    script:
        "src/diff_analysis/trait_diff_bubble.R"

rule trait_heatmap:
    # Draw heatmap for derived traits.
    input:
        FILTERED_DERIVED_TRAITS,
        GROUPS,
        "results/data/diff_analysis/trait_post_hoc.csv"
    output:
        "results/figures/diff_analysis/trait_heatmap.pdf"
    script:
        "src/diff_analysis/trait_heatmap.R"
    

# ==================== Glycan Coexpression Module ====================
rule glycan_coexpression_clustering:
    # Perform clustering to get glycan coexpression modules.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS,
        "results/data/diff_analysis/glycan_ancova.csv"
    output:
        directory("results/figures/glycan_coexpr/cc_result/"),
        "results/data/glycan_coexpr/glycan_clusters.csv"
    script:
        "src/glycan_coexpr/cluster_glycans.R"

rule glycan_cluster_heatmap:
    # Draw heatmap for differential glycans.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        clusters="results/data/glycan_coexpr/glycan_clusters.csv"
    output:
        "results/figures/glycan_coexpr/cluster_glycan_heatmap.pdf"
    script:
        "src/glycan_coexpr/cluster_heatmap.R"

rule glycan_cor_within_cluster:
    # Correlation analysis for glycans within the same GCM or with different GCMs.
    input:
        PROCESSED_ABUNDANCE,
        "results/data/glycan_coexpr/glycan_clusters.csv"
    output:
        "results/figures/glycan_coexpr/cor_inter_intra_GCM.pdf"
    script:
        "src/glycan_coexpr/glycan_cor_within_cluster.R"

rule eigen_glycans:
    # Calculate the eigen value for each glycan cluster.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS,
        "results/data/glycan_coexpr/glycan_clusters.csv"
    output:
        "results/data/glycan_coexpr/eigen_glycans.csv"
    script:
        "src/glycan_coexpr/eigen_glycans.R"

rule cluster_ancova:
    # ANCOVA for glycan clusters.
    input:
        "results/data/glycan_coexpr/eigen_glycans.csv",
        GROUPS,
        CLINICAL
    output:
        "results/data/glycan_coexpr/cluster_ancova.csv",
        "results/data/glycan_coexpr/cluster_post_hoc.csv"
    script:
        "src/glycan_coexpr/cluster_ancova.R"

rule plot_cluster_trends:
    # Plot the alteration trends of glycan clusters from the heatmap about.
    input:
        "results/data/glycan_coexpr/eigen_glycans.csv",
        GROUPS,
        "results/data/glycan_coexpr/cluster_post_hoc.csv"
    output:
        "results/figures/glycan_coexpr/glycan_cluster_trends.pdf"
    script:
        "src/glycan_coexpr/plot_cluster_trends.R"

rule plot_cluster_properties:
    # Plot heatmap showing glycans' properties.
    input:
        "results/data/glycan_coexpr/glycan_clusters.csv",
        META_PROPERTIES
    output:
        "results/figures/glycan_coexpr/glycan_property_heatmap.pdf"
    script:
        "src/glycan_coexpr/cluster_property_heatmap.R"

rule cluster_corrplot:
    # Draw cluster corrplot.
    input:
        "results/data/glycan_coexpr/eigen_glycans.csv"
    output:
        "results/figures/glycan_coexpr/cluster_corrplot.pdf"
    script:
        "src/glycan_coexpr/cluster_cor.R"

# ==================== Correlation with Liver Function ====================
rule clinical_subtype_boxplots:
    # Draw boxplots for differential traits in clinical subtypes.
    input:
        FILTERED_DERIVED_TRAITS,
        GROUPS,
        CLINICAL
    output:
        "results/figures/cor_with_liver_function/trait_clinical_subtype_boxplots.pdf"
    script:
        "src/cor_with_liver_function/clinical_subtype_boxplots.R"

rule relation_with_liver_function:
    # Perform statistical analysis to reveal relations between glycosylation and
    # liver function-related clinical variables.
    input:
        PROCESSED_ABUNDANCE,
        FILTERED_DERIVED_TRAITS,
        "results/data/glycan_coexpr/eigen_glycans.csv",
        GROUPS,
        CLINICAL
    output:
        "results/data/cor_with_liver_function/global_cor_result_with_liver_functions.csv",
        "results/data/cor_with_liver_function/grouped_cor_result_with_liver_functions.csv",
        "results/data/cor_with_liver_function/global_ttest_result_with_liver_functions.csv",
        "results/data/cor_with_liver_function/grouped_ttest_result_with_liver_functions.csv"
    script:
        "src/cor_with_liver_function/relation_with_liver_function.R"

rule liver_function_heatmap:
    # Draw the heatmap showing the relationship between glycans and
    # liver function clinical variables.
    input:
        "results/data/cor_with_liver_function/global_cor_result_with_liver_functions.csv",
        "results/data/cor_with_liver_function/global_ttest_result_with_liver_functions.csv"
    output:
        "results/figures/cor_with_liver_function/glycan_liver_function_heatmap.pdf",
        "results/figures/cor_with_liver_function/trait_liver_function_heatmap.pdf"
    script:
        "src/cor_with_liver_function/glycan_with_clinical_heatmap.R"

rule glycan_liver_function_cor_coef_per_group_boxplot:
    # Draw boxplots to show Corr Coef difference for four groups.
    input:
        "results/data/cor_with_liver_function/grouped_cor_result_with_liver_functions.csv"
    output:
        "results/figures/cor_with_liver_function/cor_coef_per_group.pdf"
    script:
        "src/cor_with_liver_function/cor_coef_per_group_boxplot.R"
        
rule ALBI_model:
    # Train a model using glycans to predict ALBI stages.
    input:
        PROCESSED_ABUNDANCE,
        CLINICAL,
        GROUPS
    output:
        "results/data/cor_with_liver_function/ALBI_model_scores.csv",
        "results/data/cor_with_liver_function/ALBI_model_preds.csv"
    script:
        "src/cor_with_liver_function/ALBI_model.R"

rule ALBI_model_ROC:
    # Draw ROC curve for ALBI model.
    input:
        "results/data/cor_with_liver_function/ALBI_model_preds.csv",
        "results/data/cor_with_liver_function/ALBI_model_scores.csv"
    output:
        "results/figures/cor_with_liver_function/ALBI_model_ROC.pdf"
    script:
        "src/cor_with_liver_function/ALBI_model_ROC.R"

rule ALBI_model_confusion_matrix:
    # Draw confusion matrices for ALBI model.
    input:
        "results/data/cor_with_liver_function/ALBI_model_preds.csv",
        GROUPS
    output:
        "results/figures/cor_with_liver_function/ALBI_model_confusion_matrix.pdf",
        "results/figures/cor_with_liver_function/ALBI_model_confusion_matrix_grouped.pdf"
    script:
        "src/cor_with_liver_function/ALBI_model_confusion_matrix.R"

rule ALBI_model_prob_ridge_plot:
    # Draw probability ridge plot for ALBI model.
    input:
        "results/data/cor_with_liver_function/ALBI_model_preds.csv",
        GROUPS
    output:
        "results/figures/cor_with_liver_function/ALBI_model_prob_ridge_plot.pdf"
    script:
        "src/cor_with_liver_function/ALBI_model_prob_ridge_plot.R"

rule ALBI_score_model:
    # Train a regression model to predict the ALBI score.
    input:
        PROCESSED_ABUNDANCE,
        CLINICAL,
        GROUPS
    output:
        "results/data/cor_with_liver_function/ALBI_score_model_metrics.csv",
        "results/figures/cor_with_liver_function/ALBI_score_model_prediction.pdf"
    script:
        "src/cor_with_liver_function/ALBI_score_model.R"


# ==================== Molecular Subtypes ====================
rule consensus_clustering:
    # Perform consensus clustering on glycans.
    input:
        RAW_ABUNDANCE,
        PROCESSED_ABUNDANCE,
        GROUPS,
        "data/plates.csv"
    output:
        directory("results/figures/subtypes/cc_result/"),
        "results/data/subtypes/consensus_cluster_result.csv",
    script:
        "src/subtypes/consensus_clustering.R"

rule subtype_pca:
    # Draw PCA plot for subtypes.
    input:
        PROCESSED_ABUNDANCE,
        "results/data/subtypes/consensus_cluster_result.csv"
    output:
        "results/figures/subtypes/subtype_pca.pdf"
    script:
        "src/subtypes/subtype_pca.R"

rule subtype_glycan_diff:
    # Perform differential analysis on glycans between subtypes.
    input:
        PROCESSED_ABUNDANCE,
        "results/data/subtypes/consensus_cluster_result.csv"
    output:
        "results/data/subtypes/subtype_glycan_anova.csv",
        "results/data/subtypes/subtype_glycan_post_hoc.csv"
    script:
        "src/subtypes/subtype_glycan_diff.R"

rule subtype_heatmap:
    # Draw heatmap for glycans in different subtypes.
    input:
        abundance=PROCESSED_ABUNDANCE,
        clusters="results/data/subtypes/consensus_cluster_result.csv",
        anova_result="results/data/subtypes/subtype_glycan_anova.csv",
        coexp_modules="results/data/glycan_coexpr/glycan_clusters.csv",
        clinical=CLINICAL
    output:
        "results/figures/subtypes/subtype_heatmap.pdf"
    script:
        "src/subtypes/subtype_heatmap.R"

rule subtype_with_clinical:
    # Inquire the relationship between glycan molecular subtypes with
    # clinical variables including liver function markers and
    # clinical stagings.
    input:
        "results/data/subtypes/consensus_cluster_result.csv",
        CLINICAL
    output:
        "results/data/subtypes/subtype_with_continous_clinical_kruskal_result.csv",
        "results/data/subtypes/subtype_with_continous_clinical_post_hoc_result.csv",
        "results/data/subtypes/subtype_with_categoric_clinical_fisher_result.csv",
        "results/data/subtypes/subtype_with_categoric_clinical_post_hoc_result.csv"
    script:
        "src/subtypes/subtype_with_clinical.R"

rule subtype_coexp_module_diff:
    # Perform differential analysis on glycan coexpression modules between subtypes.
    input:
        "results/data/glycan_coexpr/eigen_glycans.csv",
        "results/data/subtypes/consensus_cluster_result.csv"
    output:
        "results/data/subtypes/subtype_coexp_module_anova.csv",
        "results/data/subtypes/subtype_coexp_module_post_hoc.csv"
    script:
        "src/subtypes/subtype_coexp_module_diff.R"

rule subtype_coexp_module_heatmap:
    # Draw heatmap for glycan coexpression modules in different subtypes.
    input:
        "results/data/glycan_coexpr/eigen_glycans.csv",
        "results/data/subtypes/consensus_cluster_result.csv",
    output:
        "results/figures/subtypes/subtype_coexp_module_heatmap.pdf"
    script:
        "src/subtypes/subtype_coexp_module_heatmap.R"

rule subtype_continuous_clinical_boxplot:
    # Draw boxplots for continuous clinical variables in different subtypes.
    input:
        CLINICAL,
        "results/data/subtypes/consensus_cluster_result.csv",
        "results/data/subtypes/subtype_with_continous_clinical_post_hoc_result.csv"
    output:
        "results/figures/subtypes/subtype_continuous_clinical_boxplots.pdf"
    script:
        "src/subtypes/subtype_continuous_clinical_boxplot.R"

rule hypergeometric_test_heatmap:
    # Draw heatmap for hypergeometric test results.
    input:
        "results/data/subtypes/subtype_with_categoric_clinical_fisher_result.csv",
        "results/data/subtypes/subtype_with_categoric_clinical_post_hoc_result.csv"
    output:
        "results/figures/subtypes/hypergeometric_test_heatmap.pdf"
    script:
        "src/subtypes/hypergeometric_test_heatmap.R"

rule fisher_test_barplot:
    # Draw barplot for Fisher test results.
    input:
        clinical=CLINICAL,
        subtypes="results/data/subtypes/consensus_cluster_result.csv",
        fisher_result="results/data/subtypes/subtype_with_categoric_clinical_fisher_result.csv"
    output:
        "results/figures/subtypes/fisher_test_barplot.pdf"
    script:
        "src/subtypes/fisher_test_barplot.R"

rule compare_with_other_group_heatmap:
    # Draw heatmap for comparing glycan subtypes with other groups.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        subtypes="results/data/subtypes/consensus_cluster_result.csv"
    output:
        "results/figures/subtypes/compare_with_other_group_heatmap.pdf"
    script:
        "src/subtypes/compare_with_other_group_heatmap.R"

rule pairwise_comparison_of_glycans_between_subtypes_and_other_groups:
    # Perform pairwise t-test between subtypes and other groups
    # to find differential glycans.
    # Also, draw heatmap for comparing glycan subtypes with other groups.
    # This heatmap shows the fold changes of each comparison,
    # as well as the significance of the difference.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        subtypes="results/data/subtypes/consensus_cluster_result.csv"
    output:
        "results/data/subtypes/subtype_with_other_groups_glycan_ttest.csv",
        "results/figures/subtypes/compare_with_other_group_glycan_diff_heatmap.pdf"
    script:
        "src/subtypes/subtype_control_pairwise_glycan_diff.R"

rule pairwise_comparison_of_clinical_between_subtypes_and_other_groups:
    # Perform pairwise t-test between subtypes and other groups
    # to find differential clinical variables.
    input:
        clinical=CLINICAL,
        groups=GROUPS,
        subtypes="results/data/subtypes/consensus_cluster_result.csv"
    output:
        "results/data/subtypes/subtype_with_other_groups_clinical_wilcox.csv",
        "results/figures/subtypes/compare_with_other_group_clinical_diff_heatmap.pdf"
    script:
        "src/subtypes/subtype_control_pairwise_clinical_diff.R"


# ==================== TCGA Gene Expression ====================
rule download_TCGA:
    # Download TCGA gene expression data.
    output:
        "results/data/TCGA/query_summary.csv",
        directory("results/data/GDCdata"),
        "results/data/TCGA/prepared_data.rda"
    script:
        "src/TCGA/download.R"
        
rule DEA_TCGA:
    # Perform differential expression analysis on TCGA data.
    input:
        "results/data/TCGA/prepared_data.rda"
    output:
        "results/data/TCGA/dea_results.csv"
    script:
        "src/TCGA/DEA.R"

rule volcano_TCGA:
    # Draw volcano plot for TCGA gene expression.
    input:
        "results/data/TCGA/dea_results.csv",
        "data/glycogenes.csv"
    output:
        "results/figures/TCGA/volcano.pdf"
    script:
        "src/TCGA/volcano.R"

rule MA_TCGA:
    # Draw MA plot for TCGA gene expression.
    input:
        "results/data/TCGA/dea_results.csv",
        "data/glycogenes.csv"
    output:
        "results/figures/TCGA/MA_plot.pdf"
    script:
        "src/TCGA/MA_plot.R"

rule heatmap_TCGA:
    # Draw heatmap for differential glycogenes.
    input:
        "results/data/TCGA/prepared_data.rda",
        "results/data/TCGA/dea_results.csv",
        "data/glycogenes.csv"
    output:
        "results/figures/TCGA/heatmap.pdf"
    script:
        "src/TCGA/heatmap.R"

rule TCGA_consensus_cluster:
    # Perform consensus clustering using glycogenes on tumor samples.
    input:
        "results/data/TCGA/prepared_data.rda",
        "data/glycogenes.csv"
    output:
        directory("results/figures/TCGA/consensus_cluster/"),
        "results/data/TCGA/consensus_cluster_result.csv"
    script:
        "src/TCGA/consensus_cluster.R"

rule DEA_cluster:
    # DEA on the two clusters identified by consensus clustering.
    input:
        "results/data/TCGA/prepared_data.rda",
        "results/data/TCGA/consensus_cluster_result.csv"
    output:
        "results/data/TCGA/cluster_dea_results.csv"
    script:
        "src/TCGA/DEA_cluster.R"

rule TCGA_cluster_heatmap:
    # Draw heatmap for different clusters.
    input:
        "results/data/TCGA/prepared_data.rda",
        "results/data/TCGA/consensus_cluster_result.csv",
        "data/glycogenes.csv"
    output:
        "results/figures/TCGA/cluster_dea_heatmap.pdf"
    script:
        "src/TCGA/heatmap_cluster.R"

rule TCGA_clinical_heatmap:
    # Draw clinical heatmap for TCGA samples.
    input:
        "results/data/TCGA/prepared_data.rda"
    output:
        "results/figures/TCGA/clinical_heatmap.pdf"
    script:
        "src/TCGA/clinical_heatmap.R"

rule TCGA_single_gene_survival:
    # Perform Cox PH on each glycogenes.
    input:
        "results/data/TCGA/prepared_data.rda",
        "data/glycogenes.csv"
    output:
        "results/data/TCGA/cox.csv",
        directory("results/figures/TCGA/single_gene_KM/")
    script:
        "src/TCGA/single_gene_survival.R"

rule TCGA_cluster_KM:
    # Draw KM Curve for two clusters.
    input:
        "results/data/TCGA/prepared_data.rda",
        "results/data/TCGA/consensus_cluster_result.csv"
    output:
        "results/figures/TCGA/cluster_KM.pdf"
    script:
        "src/TCGA/cluster_survival.R"


# ==================== Machine Learning ====================


# ==================== Others ====================
rule draw_complete_glycans:
    # Draw glycans with glycowork with linkage information.
    input:
        "data/glycan_strucutre_best_resolution.csv"
    output:
        directory("results/figures/SNFG/complete_svg/")
    run:
        import csv
        from pathlib import Path
        from glycowork.motif.draw import GlycoDraw

        output_dir = output[0]
        Path(output_dir).mkdir(exist_ok=True, parents=True)
        with open(input[0], encoding='utf-8-sig') as fp:
            reader = csv.reader(fp)
            next(reader)
            for composition, structure in reader:
                path = Path(output_dir) / f"{composition}.svg"
                GlycoDraw(structure, filepath=str(path))

rule draw_compact_glycans:
    # Draw glycans with glycowork in the compact form.
    input:
        "data/glycan_structure_guess_linkage.csv"
    output:
        directory("results/figures/SNFG/compact_svg/")
    run:
        import csv
        from pathlib import Path
        from glycowork.motif.draw import GlycoDraw

        output_dir = output[0]
        Path(output_dir).mkdir(exist_ok=True, parents=True)
        with open(input[0], encoding='utf-8-sig') as fp:
            reader = csv.reader(fp)
            next(reader)
            for composition, structure in reader:
                path = Path(output_dir) / f"{composition}.svg"
                GlycoDraw(structure, filepath=str(path), compact=True)

rule convert_SNFG_svg_to_pdf:
    input:
        "results/figures/SNFG/{type}_svg/"
    output:
        directory("results/figures/SNFG/{type}_pdf/")
    run:
        from pathlib import Path
        from svglib.svglib import svg2rlg
        from reportlab.graphics import renderPDF
        output_path = Path(output[0])
        output_path.mkdir(parents=True)
        for svg_file in Path(input[0]).glob("*.svg"):
            name = svg_file.stem
            pdf_file = output_path / f"{name}.pdf"
            drawing = svg2rlg(str(svg_file))
            renderPDF.drawToFile(drawing, str(pdf_file))
