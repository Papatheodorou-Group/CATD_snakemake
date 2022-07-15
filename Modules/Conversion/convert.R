## Format conversion script for anndata -> Seurat object conversions using sceasy.
##
## @zgr2788



filename <- sub("Input/", "", commandArgs(trailingOnly = TRUE))

library(devtools)
devtools::install_github('cellgeni/sceasy')
library(sceasy)
library(Seurat)
sceasy::convertFormat(paste0("Input/", filename), from = "anndata", to = "seurat", outFile = paste0("Scratch/Conversion/", sub(".h5ad", "", filename), ".rds"))
