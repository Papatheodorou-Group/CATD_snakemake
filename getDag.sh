#!/bin/bash

source /nfs/research/irene/ozgurb/soft/mambaforge/etc/profile.d/conda.sh
conda activate snakemake
snakemake -n --dag | dot -Tpng -Gdpi=300 > dag.png
