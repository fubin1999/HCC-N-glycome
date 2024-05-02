rule all:
    input:
        "results/data/raw_abundance.csv"

rule build_db:
    input:
        "data/human_serum_glycans.csv"
    output:
        "results/data/glycan_db.byonic"
    script:
        "src/preprocess/make_db.py"

rule run_glyhunter:
    input:
        mass_list="data/mass_list/plate{no}.xlsx",
        db="results/data/glycan_db.byonic"
    output:
        directory("results/data/glyhunter_results/plate{no}/")
    shell:
        "glyhunter run {input.mass_list} -d {input.db} -o {output}"

rule assign_maldi_pos:
    input:
        "data/MALDI_plate_position/plate{no}.csv",
        "results/data/glyhunter_results/plate{no}/"
    output:
        "results/data/data_per_plate/plate{no}.csv"
    script:
        "src/preprocess/assign_maldi_pos.R"

rule combine_plates:
    input:
        expand("results/data/data_per_plate/plate{no}.csv", no=range(1, 9))
    output:
        "results/data/raw_abundance.csv"
    run:
        import pandas as pd
        dfs = [pd.read_csv(f) for f in input]
        combined = pd.concat(dfs)
        combined.to_csv(output[0], index=False)
