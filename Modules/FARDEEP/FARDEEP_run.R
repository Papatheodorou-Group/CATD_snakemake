## FARDEEP deconv script
##
## @zgr2788


#Load FARDEEP
if (!require("FARDEEP", quietly = TRUE))
    install.packages("FARDEEP", repos='http://cran.us.r-project.org')
suppressMessages(library(FARDEEP))
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
#Get results and reorder the matrices for correspondence
message('Started running')
res <- fardeep(C1, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)
res <- t(res$abs.beta)
res <- apply(res,2,function(x) x/sum(x)) #explicit STO constraint
res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
