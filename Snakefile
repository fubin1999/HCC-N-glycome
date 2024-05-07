PREPARED_DIR = "results/data/prepared/"
RAW_ABUNDANCE = PREPARED_DIR + "raw_abundance.csv"
PROCESSED_ABUNDANCE = PREPARED_DIR + "processed_abundance.csv"
GROUPS = PREPARED_DIR + "groups.csv"
CLINICAL = PREPARED_DIR + "clinical.csv"
DERIVED_TRAITS = PREPARED_DIR + "derived_traits.csv"

rule all:
    input:
        # ===== Prepared Data =====
        PROCESSED_ABUNDANCE,
        RAW_ABUNDANCE,
        GROUPS,
        CLINICAL,
        DERIVED_TRAITS,

        # ===== Data Quality Figures =====
        "results/figures/data_quality/batch_effect_pca.pdf"


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

rule derived_traits:
    # Calculate derived traits using GlyTrait
    input:
        PROCESSED_ABUNDANCE,
        "data/human_serum_glycans.csv"
    output:
        DERIVED_TRAITS
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