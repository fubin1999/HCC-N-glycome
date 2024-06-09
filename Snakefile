PREPARED_DIR = "results/data/prepared/"
RAW_ABUNDANCE = PREPARED_DIR + "raw_abundance.csv"
PROCESSED_ABUNDANCE = PREPARED_DIR + "processed_abundance.csv"
GROUPS = PREPARED_DIR + "groups.csv"
CLINICAL = PREPARED_DIR + "clinical.csv"

rule all:
    input:
        # ===== Others =====
        "results/data/clinical/AFP_cutoff.csv",
        "results/figures/clinical/AFP_cutoff.pdf",
        "results/figures/clinical/clinical_heatmap.pdf",
        "results/figures/SNFG/complete_svg/",
        "results/figures/SNFG/compact_svg/",
        "results/figures/SNFG/complete_pdf/",
        "results/figures/SNFG/compact_pdf/",

        # ===== Data Quality Figures =====
        "results/figures/data_quality/batch_effect_pca.pdf",

        # ===== Glycan Abundance Data =====
        "results/data/glycan_abundance/ancova_for_glycans.csv",
        "results/data/glycan_abundance/fold_change.csv",
        "results/data/glycan_abundance/roc_auc.csv",
        "results/data/glycan_abundance/cor_with_AFP.csv",

        # ===== Glycan Abundance Figures =====
        "results/figures/glycan_abundance/diff_rose_plot.pdf",
        "results/figures/glycan_abundance/diff_glycan_heatmap.pdf",
        "results/figures/glycan_abundance/glycan_cluster_trends.pdf",
        "results/figures/glycan_abundance/diff_bubble.pdf",
        "results/figures/glycan_abundance/violin_plots.pdf",
        "results/figures/glycan_abundance/diff_upset.pdf",
        "results/figures/glycan_abundance/confounders.pdf",
        "results/figures/glycan_abundance/pca.pdf",
        "results/figures/glycan_abundance/cor_with_AFP_1.pdf",
        "results/figures/glycan_abundance/cor_with_AFP_2.pdf",

        # ===== Derived Traits Data =====
        "results/data/derived_traits/derived_traits.csv",
        "results/data/derived_traits/meta_properties.csv",
        "results/data/derived_traits/ancova_for_derived_traits.csv",
        "results/data/derived_traits/posthoc_for_derived_traits.csv",
        "results/data/derived_traits/fold_change.csv",
        "results/data/derived_traits/AFP_subtype_ancova.csv",
        "results/data/derived_traits/corr_with_clinical.csv",

        # ===== Derived Traits Figures =====
        "results/figures/derived_traits/heatmap.pdf",
        "results/figures/derived_traits/boxplots_for_selected_traits.pdf",
        "results/figures/derived_traits/diff_antenna_trait_radar.pdf",
        "results/figures/derived_traits/diff_bubble.pdf",
        "results/figures/derived_traits/confounders.pdf",
        "results/figures/derived_traits/CFc_AFP_subtype_boxplot.pdf",
        "results/figures/derived_traits/corr_with_clinical.pdf",
        "results/figures/derived_traits/scatter_with_clinical.pdf",

        # ===== Residues Data =====
        "results/data/residues/glycan_residues.csv",
        "results/data/residues/ancova_result.csv",
        "results/data/residues/post_hoc_result.csv",

        # ===== Residues Figures =====
        "results/figures/residues/residue_heatmap.pdf",
        "results/figures/residues/residue_boxplots.pdf",

        # ===== GlyCompare Data =====
        "results/data/GlyCompare_results/",

        # ===== GlyCompare Figures =====
        "results/figures/SNFG/glycomotifs/",

        # ===== TCGA Data =====
        "results/data/TCGA/dea_results.csv",
        "results/data/TCGA/consensus_cluster_result.csv",
        "results/data/TCGA/cluster_dea_results.csv",

        # ===== TCGA Figures =====
        "results/figures/TCGA/volcano.pdf",
        "results/figures/TCGA/MA_plot.pdf",
        "results/figures/TCGA/heatmap.pdf",
        "results/figures/TCGA/consensus_cluster/",
        "results/figures/TCGA/cluster_dea_heatmap.pdf",
        "results/figures/TCGA/clinical_heatmap.pdf",
        "results/figures/TCGA/survival/",

        # ===== Machine Learning Data =====
        "results/data/ml/model_comparison.csv",
        "results/data/ml/all_combination_feature_selection_result.csv",
        "results/data/ml/predictions.csv",
        "results/data/ml/model_performance.csv",

        # ===== Machine Learning Figures =====
        "results/figures/ml/model_comparison_heatmap.pdf",
        "results/figures/ml/roc_curves.pdf",
        "results/figures/ml/pr_curves.pdf",
        "results/figures/ml/calibration_curve.pdf",
        "results/figures/ml/confusion_matrix.pdf",
        "results/figures/ml/probability_boxplots.pdf",
        "results/figures/ml/complex_model_metrics_table.pdf",
        "results/figures/ml/simple_model_metrics_table.pdf",
        "results/figures/ml/forest_plot.pdf",
        "results/figures/ml/shap_summary.pdf",
        "results/figures/ml/shap_waterfall/",
        "results/figures/ml/decision_tree.svg"


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

rule assign_maldi_pos:
    # Match the MALDI positions of GlyHunter results to samples.
    input:
        "data/MALDI_positions.csv",
        "results/data/glyhunter_results/plate{no}/"
    output:
        "results/data/data_per_plate/plate{no}.csv"
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

rule preprocess:
    # Filter glycan, impute missing values, and normalize.
    input:
        RAW_ABUNDANCE
    output:
        PROCESSED_ABUNDANCE
    script:
        "src/prepare_data/preprocess.R"

rule prepare_groups:
    # Prepare the groups.
    input:
        plates="data/plates.csv",
        abundance=PROCESSED_ABUNDANCE
    output:
        GROUPS
    script:
        "src/prepare_data/prepare_groups.R"

rule prepare_clinical:
    # Prepare the clinical information.
    input:
        clinical="data/clinical.csv",
        abundance=PROCESSED_ABUNDANCE,
        plates="data/plates.csv"
    output:
        CLINICAL
    script:
        "src/prepare_data/prepare_clinical.R"


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


# ==================== Glycan Abundance ====================
rule ancova_for_glycans:
    # Perform ANCOVA for each glycan.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        clinical=CLINICAL
    output:
        ancova="results/data/glycan_abundance/ancova_for_glycans.csv",
        posthoc="results/data/glycan_abundance/posthoc_for_glycans.csv"
    script:
        "src/glycan_abundance/ancova_for_glycans.R"

rule fold_change:
    # Calculate fold changes for each glycan.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS
    output:
        "results/data/glycan_abundance/fold_change.csv"
    script:
        "src/glycan_abundance/fold_change.R"

rule diff_rose_plot:
    # Draw rose plot for differential glycans between each group pair.
    input:
        "results/data/glycan_abundance/posthoc_for_glycans.csv"
    output:
        "results/figures/glycan_abundance/diff_rose_plot.pdf"
    script:
        "src/glycan_abundance/rose_plot.R"

rule diff_glycan_heatmap:
    # Draw heatmap for differential glycans.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        ancova_result="results/data/glycan_abundance/ancova_for_glycans.csv",
        mp_table="results/data/derived_traits/meta_properties.csv"
    output:
        "results/figures/glycan_abundance/diff_glycan_heatmap.pdf",
        "results/data/glycan_abundance/glycan_clusters.csv"
    script:
        "src/glycan_abundance/heatmap.R"

rule diff_bubble:
    # Draw bubble plot for p-values and fold changes of glycans.
    input:
        post_hoc="results/data/glycan_abundance/posthoc_for_glycans.csv",
        fold_change="results/data/glycan_abundance/fold_change.csv",
        row_order="results/data/glycan_abundance/glycan_clusters.csv"
    output:
        "results/figures/glycan_abundance/diff_bubble.pdf"
    script:
        "src/glycan_abundance/diff_bubble.R"

rule glycan_cluster_trends:
    # Plot the alteration trends of glycan clusters from the heatmap about.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        clusters="results/data/glycan_abundance/glycan_clusters.csv"
    output:
        "results/figures/glycan_abundance/glycan_cluster_trends.pdf"
    script:
        "src/glycan_abundance/cluster_trends.R"

rule violin_plots:
    # Plot violin plots for all significant glycans.
    input:
        abundance=PROCESSED_ABUNDANCE,
        groups=GROUPS,
        ancova_result="results/data/glycan_abundance/ancova_for_glycans.csv"
    output:
        "results/figures/glycan_abundance/violin_plots.pdf"
    script:
        "src/glycan_abundance/violin_plots.R"

rule diff_upset:
    # Plot upset plot for the number of significant glycans
    # between each group paires.
    input:
        "results/data/glycan_abundance/posthoc_for_glycans.csv"
    output:
        "results/figures/glycan_abundance/diff_upset.pdf"
    script:
        "src/glycan_abundance/diff_upset.R"

rule confounders:
    # Plot dot plot for confounders' p-values.
    input:
        "results/data/glycan_abundance/ancova_for_glycans.csv"
    output:
        "results/figures/glycan_abundance/confounders.pdf"
    script:
        "src/glycan_abundance/confounder_dot_plot.R"

rule glycan_roc:
    # Perform ROC analysis on glycans, and draw ROC curves.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS
    output:
        "results/data/glycan_abundance/roc_auc.csv",
        "results/figures/glycan_abundance/roc_curves.pdf"
    script:
        "src/glycan_abundance/roc.R"

rule glycan_pca:
    # Draw PCA plots for glycan abundance.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS
    output:
        "results/figures/glycan_abundance/pca.pdf"
    script:
        "src/glycan_abundance/pca.R"

rule glycan_cor_with_AFP:
    # Draw corrplot for glycans and AFP.
    input:
        PROCESSED_ABUNDANCE,
        GROUPS,
        CLINICAL
    output:
        "results/data/glycan_abundance/cor_with_AFP.csv",
        # There are too many glycans to display in a row,
        # and corrplot object could not be jointed by cowplot,
        # so the plot is splited into two.
        "results/figures/glycan_abundance/cor_with_AFP_1.pdf",
        "results/figures/glycan_abundance/cor_with_AFP_2.pdf"
    script:
        "src/glycan_abundance/corr_with_AFP.R"


# ==================== Derived Traits ====================
rule calculate_derived_traits:
    # Calculate derived traits using GlyTrait
    input:
        PROCESSED_ABUNDANCE,
        "data/human_serum_glycans.csv",
        "src/derived_traits/struc_builtin_formulas.txt"
    output:
        filtered_traits="results/data/derived_traits/filtered_derived_traits.csv",
        all_traits="results/data/derived_traits/derived_traits.csv",
        mp_table="results/data/derived_traits/meta_properties.csv"
    script:
        "src/derived_traits/derive_traits.py"

rule trait_ancova:
    # Perform ANCOVA on derived traits.
    input:
        traits="results/data/derived_traits/filtered_derived_traits.csv",
        groups=GROUPS,
        clinical=CLINICAL
    output:
        "results/data/derived_traits/ancova_for_derived_traits.csv",
        "results/data/derived_traits/posthoc_for_derived_traits.csv"
    script:
        "src/derived_traits/ancova_for_derived_traits.R"

rule trait_fold_change:
    # Calculate fold changes for each derived trait.
    input:
        "results/data/derived_traits/filtered_derived_traits.csv",
        GROUPS
    output:
        "results/data/derived_traits/fold_change.csv"
    script:
        "src/derived_traits/fold_change.R"

rule trait_confounders:
    # Plot dot plot for confounders' p-values.
    input:
        "results/data/derived_traits/ancova_for_derived_traits.csv"
    output:
        "results/figures/derived_traits/confounders.pdf"
    script:
        "src/derived_traits/confounder_dot_plot.R"

rule boxplots_for_selected_traits:
    # Draw boxplots for selected derived traits.
    input:
        "results/data/derived_traits/derived_traits.csv",
        "results/data/prepared/groups.csv"
    output:
        "results/figures/derived_traits/boxplots_for_selected_traits.pdf"
    script:
        "src/derived_traits/selected_boxplots.R"

rule diff_antenna_trait_radar:
    # Draw radar plot for differential antenna traits.
    input:
        "results/data/derived_traits/derived_traits.csv",
        "results/data/prepared/groups.csv"
    output:
        "results/figures/derived_traits/diff_antenna_trait_radar.pdf"
    script:
        "src/derived_traits/diff_antenna_radar.R"

rule trait_diff_bubble:
    # Draw bubble plot for differential derived traits.
    input:
        "results/data/derived_traits/posthoc_for_derived_traits.csv"
    output:
        "results/figures/derived_traits/diff_bubble.pdf"
    script:
        "src/derived_traits/diff_bubble.R"

rule AFP_subtype_trait_diff:
    # Perform ANCOVA on derived traits between AFP negative and AFP positive HCC samples.
    input:
        "results/data/derived_traits/filtered_derived_traits.csv",
        GROUPS,
        CLINICAL
    output:
        "results/data/derived_traits/AFP_subtype_ancova.csv",
        "results/figures/derived_traits/CFc_AFP_subtype_boxplot.pdf"
    script:
        "src/derived_traits/AFP_subtype_ancova.R"

rule trait_heatmap:
    # Draw heatmap for derived traits.
    input:
        "results/data/derived_traits/filtered_derived_traits.csv",
        GROUPS,
        "results/data/derived_traits/posthoc_for_derived_traits.csv"
    output:
        "results/figures/derived_traits/heatmap.pdf"
    script:
        "src/derived_traits/heatmap.R"

rule trait_cor_with_clinical:
    # Draw corrplot for derived traits and clinical data.
    input:
        "results/data/derived_traits/filtered_derived_traits.csv",
        GROUPS,
        CLINICAL
    output:
        "results/data/derived_traits/corr_with_clinical.csv",
        "results/data/derived_traits/corr_with_clinical_HCC.csv",
        "results/figures/derived_traits/corr_with_clinical.pdf",
        "results/figures/derived_traits/corr_with_clinical_HCC.pdf"
    script:
        "src/derived_traits/corr_with_clinical.R"

rule scatter_with_clinical:
    # Draw scatter plot for selected derived traits and clinical variables
    # with high correlation.
    input:
        "results/data/derived_traits/filtered_derived_traits.csv",
        GROUPS,
        CLINICAL
    output:
        "results/figures/derived_traits/scatter_with_clinical.pdf"
    script:
        "src/derived_traits/scatter_with_clinical.R"


# ==================== Residue Analysis ====================
rule calculate_residues:
    # Calculate glycan residues (the abundance of each monosaccharide).
    input:
        PROCESSED_ABUNDANCE,
        "results/data/derived_traits/meta_properties.csv"
    output:
        "results/data/residues/glycan_residues.csv"
    script:
        "src/residues/calculate_residues.R"

rule ancova_for_residues:
    # Perform ANCOVA analysis on glycan residues.
    input:
        "results/data/residues/glycan_residues.csv",
        GROUPS,
        CLINICAL
    output:
        "results/data/residues/ancova_result.csv",
        "results/data/residues/post_hoc_result.csv"
    script:
        "src/residues/ancova.R"

rule residue_heatmap:
    # Draw heatmap for mean number of residues per group.
    input:
        "results/data/residues/glycan_residues.csv",
        GROUPS
    output:
        "results/figures/residues/residue_heatmap.pdf"
    script:
        "src/residues/mean_heatmap.R"

rule residue_boxplots:
    # Draw boxplots for residues.
    input:
        "results/data/residues/glycan_residues.csv",
        GROUPS
    output:
        "results/figures/residues/residue_boxplots.pdf"
    script:
        "src/residues/boxplot.R"


# ==================== GlyCompare ====================
GLYCOMPARE_CLI = "/Users/fubin/Python/glyCompareCT/glyCompareCT"

rule prepare_for_glycompare:
    # Prepare data format for GlyCompareCT.
    input:
        "data/glycan_structure_guess_linkage.csv"
    output:
        temp("results/data/glycompare_structures.csv")
    run:
        import pandas as pd
        df = pd.read_csv(input[0])
        df = df.rename(columns={"composition": "Name", "structure": "Glycan Structure"})
        df.to_csv(output[0], index=False)

rule run_glycompare:
    # Run GlyCompareCT.
    input:
        PROCESSED_ABUNDANCE,
        "results/data/glycompare_structures.csv"
    output:
        directory("results/data/GlyCompare_results/")
    shell:
        "{GLYCOMPARE_CLI} structure -a {input[0]} -v {input[1]} -o {output[0]} -p glycoCT -r N -c 8"

rule plot_motif_SNFG:
    # Draw SNFG cartoons for GlyCompare motifs.
    input:
        "results/data/GlyCompare/GlyCompare_output_data/GlyCompare_motif_annotation.csv"
    output:
        directory("results/figures/SNFG/glycomotifs/")
    run:
        import csv
        from pathlib import Path
        from glycowork.motif.draw import GlycoDraw

        output_dir = output[0]
        Path(output_dir).mkdir(exist_ok=True, parents=True)
        with open(input[0], encoding='utf-8-sig') as fp:
            reader = csv.reader(fp)
            next(reader)
            for motif, _, _, structure in reader:
                path = Path(output_dir) / f"{motif}.svg"
                GlycoDraw(structure, filepath=str(path), compact=True)


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

rule consensus_cluster:
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

rule TCGA_HCC_clinical_heatmap:
    # Draw clinical heatmap for TCGA HCC samples.
    input:
        "results/data/TCGA/prepared_data.rda",
        "results/data/TCGA/consensus_cluster_result.csv"
    output:
        "results/figures/TCGA/clinical_heatmap.pdf"
    script:
        "src/TCGA/clinical_heatmap_cluster.R"

rule TCGA_survival:
    # Draw survival plots for glycogenes.
    input:
        "results/data/TCGA/prepared_data.rda",
        "data/glycogenes.csv"
    output:
        directory("results/figures/TCGA/survival/")
    script:
        "src/TCGA/survival.R"


# ==================== Machine Learning ====================
rule prepare_data_for_ml:
    # Prepare data for ML and split the dataset.
    input:
        abundance=PROCESSED_ABUNDANCE,
        clinical=CLINICAL,
        groups=GROUPS
    output:
        train_data="results/data/ml/train_data.csv",
        test_data="results/data/ml/test_data.csv",
        feature_types="results/data/ml/feature_types.json"
    script:
        "src/ml/prepare_data.py"

rule compare_models:
    # Compare different machine learning models (including the HCC Fusion Classifier)
    # using 10-fold cross-validation.
    input:
        train_data="results/data/ml/train_data.csv",
        feature_types="results/data/ml/feature_types.json"
    output:
        "results/data/ml/model_comparison.csv"
    script:
        "src/ml/compare_models.py"

rule plot_compare_model_heatmap:
    # Draw the heatmap for model comparison
    input:
        "results/data/ml/model_comparison.csv"
    output:
        "results/figures/ml/model_comparison_heatmap.pdf"
    script:
        "src/ml/model_compare_heatmap.R"

rule all_combination_feature_selection:
    # Feature selection using all combinations of features
    input:
        "results/data/ml/train_data.csv",
        "results/data/glycan_abundance/glycan_clusters.csv"
    output:
        "results/data/ml/all_combination_feature_selection_result.csv"
    script:
        "src/ml/all_combination_feature_selection.py"

rule make_predictions:
    # Predict on the test data and evalute the model.
    input:
        train_data="results/data/ml/train_data.csv",
        test_data="results/data/ml/test_data.csv",
        feature_types="results/data/ml/feature_types.json"
    output:
        "results/data/ml/predictions.csv"
    script:
        "src/ml/make_predictions.py"

rule model_roc_and_pr_curves:
    # Plot ROC curves and PR curves for two models and AFP.
    input:
        predictions="results/data/ml/predictions.csv",
        groups=GROUPS,
        clinical=CLINICAL
    output:
        "results/figures/ml/roc_curves.pdf",
        "results/figures/ml/pr_curves.pdf"
    script:
        "src/ml/pr_roc.R"

rule confusion_matrix:
    # Draw confusion matrix for the model.
    input:
        "results/data/ml/predictions.csv",
        "results/data/prepared/groups.csv"
    output:
        "results/figures/ml/confusion_matrix.pdf"
    script:
        "src/ml/confusion_matrix.R"

rule calibration_curve:
    # Draw the calibration curve of the model.
    input:
        "results/data/ml/predictions.csv"
    output:
        "results/figures/ml/calibration_curve.pdf"
    script:
        "src/ml/calibration_curve.R"

rule probability_boxplots:
    # Draw boxplots displaying probability distributions for each group.
    input:
        "results/data/ml/predictions.csv",
        "results/data/prepared/groups.csv"
    output:
        "results/figures/ml/probability_boxplots.pdf"
    script:
        "src/ml/score_boxplots.R"

rule evaluate_model:
    # Calculate the performance metrics of the model.
    input:
        "results/data/ml/predictions.csv",
        "results/data/prepared/groups.csv"
    output:
        "results/data/ml/model_performance.csv"
    script:
        "src/ml/metrics.R"

rule model_metrics_table:
    # Draw a table of model metrics.
    input:
        "results/data/ml/model_performance.csv"
    output:
        "results/figures/ml/complex_model_metrics_table.pdf",
        "results/figures/ml/simple_model_metrics_table.pdf"
    script:
        "src/ml/metrics_table.R"

rule forest_plot:
    # Draw a forest plot for the HCC Slim model.
    input:
        "results/data/ml/train_data.csv"
    output:
        "results/figures/ml/forest_plot.pdf"
    script:
        "src/ml/forest_plot.R"

rule shap:
    # SHAP on HCC Slim.
    input:
        "results/data/ml/train_data.csv",
        "results/data/ml/test_data.csv"
    output:
        "results/figures/ml/shap_summary.pdf",
        directory("results/figures/ml/shap_waterfall/")
    script:
        "src/ml/shap.py"

rule decision_tree:
    # Visualize a decision tree classifier.
    input:
        "results/data/ml/train_data.csv",
        "results/data/ml/test_data.csv"
    output:
        "results/figures/ml/decision_tree.svg"
    notebook:
        "src/ml/decision_tree.ipynb"


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