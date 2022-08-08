#!/bin/bash

#BSUB -n 64
#BSUB -M 307200
#BSUB -R "rusage[mem=307200]"
#BSUB -o "pipOut.log"
#BSUB -J "runPip"

condaDir="/nfs/research/irene/ozgurb/soft/mambaforge/etc/profile.d/conda.sh"

source $condaDir
conda activate snakemake
snakemake \
	--use-conda \
	--cores 64 \
	--conda-frontend mamba
