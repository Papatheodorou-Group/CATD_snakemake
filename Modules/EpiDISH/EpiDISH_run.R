## EpiDISH deconv script
##
## @zgr2788


#Load EpiDISH
suppressMessages(library(energy))
suppressMessages(library(EpiDISH))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1 <- args[3]
filename_O <- args[4]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)


#Toss out the genes tossed out in T normalization from C as well
common <- intersect(rownames(C1), rownames(T))
C1 <- C1[common,]
T  <- T[common,]
C1 <- as.matrix(C1) #Explicit typecasting needed

#Get res and reorder the matrices for correspondence
res <- t(EpiDISH::epidish(beta.m = T, ref.m = C1, method = "RPC")$estF)
res <- apply(res,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
res <- apply(res,2,function(x) x/sum(x)) #explicit STO constraint
res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
