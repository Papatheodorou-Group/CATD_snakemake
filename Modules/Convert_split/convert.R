## Format conversion script for anndata -> Seurat object conversions using sceasy.
##
## @zgr2788



filename <- commandArgs(trailingOnly = TRUE)

suppressMessages(library(devtools))
devtools::install_github('cellgeni/sceasy',force=TRUE)
suppressMessages(library(sceasy))
use_python("/nfs/research/irene/annavp/github/CATD_snakemake/.snakemake/conda/399d2e6af93db0037d0ab81ebb046a91_/bin/python")
suppressMessages(library(Seurat))
sceasy::convertFormat(filename, from = "anndata", to = "seurat", outFile = paste0(sub(".h5ad", "", filename), "_seurat.rds"))
