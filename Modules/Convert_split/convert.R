## Format conversion script for anndata -> Seurat object conversions using sceasy.
##
## @zgr2788



filename <- commandArgs(trailingOnly = TRUE)

suppressMessages(library(devtools))
devtools::install_github('cellgeni/sceasy')
suppressMessages(library(sceasy))
suppressMessages(library(Seurat))
sceasy::convertFormat(filename, from = "anndata", to = "seurat", outFile = paste0(sub(".h5ad", "", filename), "_seurat.rds"))
