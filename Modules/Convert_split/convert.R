## Format conversion script for anndata -> Seurat object conversions using sceasy.
##
## @zgr2788



filename <- commandArgs(trailingOnly = TRUE)

library(devtools)
devtools::install_github('cellgeni/sceasy')
library(sceasy)
library(Seurat)
sceasy::convertFormat(filename, from = "anndata", to = "seurat", outFile = paste0(sub(".h5ad", "", filename), "_seurat.rds"))
