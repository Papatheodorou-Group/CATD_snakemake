#!/bin/bash

condaDir="/nfs/research/irene/ozgurb/soft/mambaforge/etc/profile.d/conda.sh"

source $condaDir
conda activate snakemake
snakemake -n --dag | dot -Tpng -Gdpi=300 > dag.png
