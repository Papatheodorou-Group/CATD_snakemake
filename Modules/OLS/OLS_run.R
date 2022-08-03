## Ordinary Least Squares deconv script
##
## @zgr2788

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

#Run OLS and reorder for correspondence
res <- apply(T,2,function(x) lm(x ~ as.matrix(C1))$coefficients[-1])
rownames(res) <- unlist(lapply(strsplit(rownames(res),")"),function(x) x[2])) #Correct rownames
res[is.na(res)] <- 0       #Convert NA's to zeros prior to applying constraints
res = apply(res,2,function(x) ifelse(x < 0, 0, x)) #Explicit non-negativity constraint
res = apply(res,2,function(x) x/sum(x)) #Explicit STO constraint
res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
