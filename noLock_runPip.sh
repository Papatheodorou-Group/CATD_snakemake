#!/bin/bash

#BSUB -o "pipOut_%J.log"
#BSUB -J "runPip"

condaDir="Run basicSetup.sh to configure"

source $condaDir
conda activate snakemake

snakemake \
  --profile lsf \
  --conda-frontend mamba \
  --nolock
