include: "Modules/Convert_split/Snakefile"
include: "Modules/Psuedobulk/Snakefile"

rule all:
    input:
        "Input/Psuedobulks/sciPlex3_MCF7_Targets45toEND_and_vehicle_pbulks.rds",
