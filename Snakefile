include: "Modules/Convert_split/Snakefile"

rule all:
    input:
        "Input/Cell_splits/sciPlex1_C0.rds",
        "Input/Cell_splits/sciPlex1_gen.rds"
