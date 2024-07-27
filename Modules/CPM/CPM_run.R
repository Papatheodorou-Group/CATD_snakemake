## CPM deconv script
##
## @zgr2788


#Load CPM
suppressMessages(library(scBio))
suppressMessages(library(energy))
library(Matrix)


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C0 <- args[3]
filename_phenData <- args[4]
cores <- as.numeric(args[5])
filename_O <- args[6]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C0 <- readRDS(filename_C0)
phenData <- readRDS(filename_phenData)
print(filename_O)

#Match genes in rows for both references
common <- intersect(rownames(C0), rownames(T))
T <- T[common,]
C0 <- C0[common,]
rm(common)

#Preprocess
#T <- T - rowMeans(T)

p <- prcomp(t(C0)[, which(apply(t(C0), 2, var) != 0)], center = TRUE, scale. = TRUE)$x[,1:2]
cellTypes <- as.character(phenData$cellType)


#Get res and reorder the matrices for correspondence
#dense matrix
C0<-as.data.frame(as.matrix(C0))


res <- t(CPM(C0, cellTypes, T, p, quantifyTypes = TRUE,typeTransformation=TRUE)$cellTypePredictions)

message("CPM running DONE")
#res = apply(res,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
#res = apply(res,2,function(x) x/sum(x)) #explicit STO constrain

if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res = res[order(match(rownames(res), rownames(P))),]

saveRDS(res, file=filename_O)
