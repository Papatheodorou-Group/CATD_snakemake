## CPM deconv script
##
## @zgr2788


#Load CPM
suppressMessages(library(scBio))
suppressMessages(library(energy))



#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C0 <- args[3]
filename_phenData <- args[4]
cores <- as.numeric(args[5])

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C0 <- readRDS(filename_C0)
phenData <- readRDS(filename_phenData)

#Match genes in rows for both references
common <- intersect(rownames(C0), rownames(T))
T <- T[common,]
C0 <- C0[common,]
rm(common)

#Preprocess
T <- T - rowMeans(T)
message("Running prcomp...")
p <- prcomp(t(C0)[, which(apply(t(C0), 2, var) != 0)], center = TRUE, scale. = TRUE)$x[,1:2]
cellTypes <- as.character(phenData$cellType)


#Get res and reorder the matrices for correspondence
res <- t(CPM(C0, cellTypes, T, p, quantifyTypes = TRUE, no_cores = cores)$cellTypePredictions)
res <- apply(res,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
res <- apply(res,2,function(x) x/sum(x)) #explicit STO constraint
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
