## dtangle deconv script
##
## @zgr2788


#Load dtangle
suppressMessages(library(dtangle))
suppressMessages(library(energy))
suppressMessages(library(tidyr))



#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1 <- args[3]
filename_C2 <- args[4]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)
C2 <- readRDS(filename_C2)


#Toss out the genes tossed out in T normalization from C as well
C1 <- C1[rownames(C1) %in% rownames(T),]
C2 <- C2[rownames(C2) %in% rownames(T),]

#Preprocess
mixtures <- t(T)
refs <- t(C1)
C2 <- tidyr::separate_rows(C2,"cellType",sep="\\|")
C2 <- C2[C2$geneName %in% rownames(C1),]
MD <- tapply(C2$geneName,C2$cellType,list)
MD <- lapply(MD,function(x) sapply(x, function(y) which(y==rownames(C1))))

#Get results and reorder the matrices for correspondence
res <- t(dtangle::dtangle(Y=mixtures, reference=refs, markers = MD)$estimates)
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
