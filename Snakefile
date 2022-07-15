include: "Modules/Convert_split/Snakefile"

rule all:
    input:
        "Output/Cell_splits/sciPlex1_C0.rds",
        "Output/Cell_splits/sciPlex1_gen.rds"
