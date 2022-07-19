configfile: 'config.yaml'
include: "Modules/Convert_split/Snakefile"
include: "Modules/Psuedobulk/Snakefile"

rule all:
    input:
        "Input/Psuedobulks/Hrvatin_afteint_pbulks.rds",
        "Input/Psuedobulks/Hrvatin_afteint_props.rds"
