#!/bin/bash

condaDir="Run basicSetup.sh to configure"

source $condaDir
conda activate snakemake
snakemake -n --dag | dot -Tpng -Gdpi=300 > dag.png
