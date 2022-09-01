## EPIC deconv script
##
## @zgr2788


#Load EPIC
suppressMessages(library(devtools))
devtools::install_github("GfellerLab/EPIC", auth_token = "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")
suppressMessages(library(EPIC))
suppressMessages(library(energy))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1 <- args[3]
filename_C2 <- args[4]
filename_refVar <- args[5]
filename_O <- args[6]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)
C2 <- readRDS(filename_C2)
refVar <- readRDS(filename_refVar)


#Toss out the genes tossed out in T normalization from C as well
common <- intersect(rownames(C1), rownames(T))
common <- intersect(common, rownames(C2))
common <- intersect(common, rownames(refVar))
C1 <- C1[common,]
C2 <- C2[common,]
T  <- T[common,]
refVar <- refVar[common,]

#Preprocess references for EPIC
markers <- as.character(C2$geneName)
cellTypes <- colnames(C1)
C_EPIC <- list()
C_EPIC[["sigGenes"]] <- rownames(C1[markers, cellTypes])
C_EPIC[["refProfiles"]] <- as.matrix(C1[markers, cellTypes])
C_EPIC[["refProfiles.var"]] <- refVar[markers, cellTypes]



#Get results and reorder the matrices for correspondence
res <- t(EPIC(bulk=as.matrix(T), reference=C_EPIC, withOtherCells=TRUE, scaleExprs=FALSE)$cellFractions) #scaleExprs=TRUE by default: only keep genes in common between matrices
res <- res[!rownames(res) %in% "otherCells",]
res[is.na(res)] <- 0
res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
