## SCDC deconv script
##
## @zgr2788


#Load SCDC
if (!require("L1pack", quietly = TRUE)){
  install.packages("Modules/SCDC/fastmatrix_0.4.tar.gz", repos=NULL, type = "source")
  install.packages("Modules/SCDC/L1pack_0.40.tar.gz", repos=NULL, type = "source")

}

suppressMessages(library(remotes))
suppressMessages(library(devtools))
remotes::install_github("renozao/xbioc")
devtools::install_github("meichendong/SCDC")
suppressMessages(library(SCDC))
suppressMessages(library(Biobase))



#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C0 <- args[3]
filename_phenData <- args[4]
filename_O <- args[5]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C0 <- readRDS(filename_C0)
phenData <- readRDS(filename_phenData)


#Match genes in rows for both references
common <- intersect(rownames(C0), rownames(T))
T <- T[common,]
C0 <- C0[common,]
rm(common)

#Convert to eset
T <- ExpressionSet(T)
C0 <- ExpressionSet(C0, phenoData = as(phenData, "AnnotatedDataFrame"))


#Get results and reorder the matrices for correspondence
res <- t(SCDC::SCDC_prop(bulk.eset = T, sc.eset = C0, ct.varname = "cellType", sample = "sampleID", ct.sub = levels(phenData$cellType), iter.max = 200)$prop.est.mvw)
res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
