## DWLS deconv script
##
## @zgr2788


#Load DWLS
suppressMessages(library(future.apply))
suppressMessages(library(DWLS))
suppressMessages(library(energy))





#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C0 <- args[3]
filename_phenData <- args[4]
plan('multisession', workers = as.numeric(args[5])) #Paralellism
filename_O <- args[6]
scratchDir <- args[7]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C0 <- readRDS(filename_C0)
phenData <- readRDS(filename_phenData)


#Match genes in rows for both references
common <- intersect(rownames(C0), rownames(T))
T <- T[common,]
C0 <- C0[common,]
rm(common)

#Preprocess
message("Started running sig...")
Signature <- buildSignatureMatrixMAST(scdata = C0, id = phenData[,"cellType"], path = scratchDir, diff.cutoff = 0.5, pval.cutoff = 0.01)


#Get results and reorder the matrices for correspondence
res <- future_apply(T,2, function(x){
  b = setNames(x, rownames(T))
  tr <- trimData(Signature, b)
  RES <- t(solveDampenedWLS(tr$sig, tr$bulk))
}, future.seed = TRUE)

rownames(res) <- as.character(unique(phenData$cellType))
res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
