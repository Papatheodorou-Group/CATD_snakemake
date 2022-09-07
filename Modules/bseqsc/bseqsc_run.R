## bseqsc deconv script
##
## @zgr2788


#Load bseqsc
suppressMessages(library(devtools))
devtools::install_github('shenorrlab/bseq-sc', auth_token = "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")
suppressMessages(library(bseqsc))
suppressMessages(library(Biobase))
suppressMessages(library(energy))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(parallel))
suppressMessages(library(preprocessCore))
suppressMessages(library(e1071))
bseqsc_config('Modules/bseqsc/CIBERSORT.R')



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
T <- data.frame(T)
C0 <- ExpressionSet(C0, phenoData = as(phenData, "AnnotatedDataFrame"))


#Preprocess
markers <- list()
for(cellType in levels(C0$cellType)){
  markers[cellType] <- list(rownames(C0))
}



#Get results and reorder the matrices for correspondence
message("Getting bseqsc_basis...")
B <- bseqsc_basis(C0, markers, clusters = 'cellType', samples = 'sampleID', ct.scale = TRUE)
message("Running CIBERSORT with B...")
res <- CIBERSORT(B, T, QN = FALSE, absolute = FALSE, perm = 0)
res <- t(res[,1:(ncol(res)-3)])
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
