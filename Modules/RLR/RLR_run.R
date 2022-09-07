## RLR deconv script
##
## @zgr2788

suppressMessages(library(energy))
suppressMessages(library(MASS))

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

#Run RLR and reorder for correspondence
message("Started running")
res <- do.call(cbind.data.frame,lapply(apply(T,2,function(x) MASS::rlm(x ~ as.matrix(C1), maxit=100)), function(y) y$coefficients[-1]))
rownames(res) <- unlist(lapply(strsplit(rownames(res),")"),function(x) x[2]))
res <- apply(res,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
res <- apply(res,2,function(x) x/sum(x)) #explicit STO constraint
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]


#Save and exit
saveRDS(res, file=filename_O)
