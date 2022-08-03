## DCQ deconv script
##
## @zgr2788


#Load ComICS

suppressMessages(library(energy))
suppressMessages(library(ComICS))


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
res <- t(dcq(reference_data = C1, mix_data = T, marker_set = as.data.frame(row.names(C1)) , alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 10)$average)
res <- apply(res,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
res <- apply(res,2,function(x) x/sum(x)) #explicit STO constraint
res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
