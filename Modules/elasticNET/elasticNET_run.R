## elasticNET deconv script
##
## @zgr2788


#Load glmnet
suppressMessages(library(glmnet))
suppressMessages(library(energy))


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
C1 <- C1[rownames(C1) %in% rownames(T),]

#Get res and reorder the matrices for correspondence
res <- apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C1), y = z, alpha = 0.2, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(C1), z)$lambda.1se))[1:ncol(C1)+1,])
res <- apply(res,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
res <- apply(res,2,function(x) x/sum(x)) #explicit STO constraint
res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
