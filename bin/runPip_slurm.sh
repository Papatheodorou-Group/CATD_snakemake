#!/bin/bash


#SBATCH --output=pipOut_%j.log
#SBATCH --job-name=runPip
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=4

condaDir=""

source $condaDir
conda activate snakemake

snakemake \
  --profile slurm \
  --conda-frontend mamba \
  --cores 32
