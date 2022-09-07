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
common <- intersect(rownames(C1), rownames(T))
C1 <- C1[common,]
T  <- T[common,]

#Run OLS and reorder for correspondence
res <- apply(T,2,function(x) lm(x ~ as.matrix(C1))$coefficients[-1])
rownames(res) <- unlist(lapply(strsplit(rownames(res),")"),function(x) x[2])) #Correct rownames
res[is.na(res)] <- 0       #Convert NA's to zeros prior to applying constraints
res = apply(res,2,function(x) ifelse(x < 0, 0, x)) #Explicit non-negativity constraint
res = apply(res,2,function(x) x/sum(x)) #Explicit STO constraint
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
