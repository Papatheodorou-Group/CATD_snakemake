## debCAM semi-supervised clustering script with cell type profile matrix (S in debCAM documentation)
##
## @zgr2788


suppressMessages(library(debCAM))
suppressMessages(library(BiocParallel))

#Read data
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_C1 <- args[2]
filename_P <- args[3]

T <- readRDS(filename_T)
C1 <- readRDS(filename_C1)
P <- readRDS(filename_P)

#Get foldchanges for possible marker genes
pMGstat <- MGstatistic(C1, colnames(C1))

#Get those with a fold change over 10, might need to revisit
pMGlist.FC <- lapply(colnames(C1), function(x) rownames(pMGstat)[pMGstat$idx == x & pMGstat$OVE.FC > 10])

#Run with marker list generated from S
res <- redoASest(T, pMGlist.FC, maxIter = 30) #Adjust maxIter later
res <- t(res$Aest)

#Convert P to matrix
P <- matrix(unlist(P), nrow = nrow(P), ncol = ncol(P))

#Calculate RMSE error
rmse <- sqrt(sum((P - res)^2) / length(res))

#Calculate euclidean distance (switch to any minkowski-type by adjusting p)
m_dist <- dist(rbind(as.vector(res), as.vector(unlist(P))), method = "minkowski", p = 2)

#print and exit (update later)
print(paste0("RMSE: ", rmse))
print(paste0("Euclidean Distance: ", m_dist))
