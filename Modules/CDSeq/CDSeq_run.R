## CDSeq deconv script
##
## @zgr2788


#Load CDseq
suppressMessages(library(devtools))
devtools::install_github("kkang7/CDSeq_R_Package")
suppressMessages(library(CDSeq))
suppressMessages(library(energy))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1 <- args[3]
cores <- as.numeric(args[4]) 

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)


#Toss out the genes tossed out in T normalization from C as well
C1 <- C1[rownames(C1) %in% rownames(T),]

#Get results and reorder the matrices for correspondence
res <- CDSeq(bulk_data = T, reference_gep = C1, cell_type_number = ncol(C1), mcmc_iterations = 1000, cpu_number = cores, block_number = 6, gene_subset_size=15)
res <- res$estProp
res <- res[order(match(rownames(res), rownames(P))),]

#Calculate RMSE error
rmse <- sqrt(mean(as.matrix((P - res)^2)))

#Calculate euclidean distance (switch to any minkowski-type by adjusting p)
m_dist <- dist(rbind(as.vector(res), as.vector(unlist(P))), method = "minkowski", p = 2)

#Calculate Distance corr
distance_corr <- dcor(P, res)

#print and exit (update later)
print(paste0("RMSE: ", rmse))
print(paste0("Euclidean Distance: ", m_dist))
print(paste0("Distance Correlation: ", distance_corr))
