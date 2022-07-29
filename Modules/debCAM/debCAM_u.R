## debCAM unsupervised deconv script
##
## @zgr2788


suppressMessages(library(debCAM))
suppressMessages(library(BiocParallel))
suppressMessages(library(energy))



#Read data
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
cellType_n <- as.numeric(args[3])


T <- readRDS(filename_T)
P <- readRDS(filename_P)


# Warning: Since fit model is linear, needs more samples in pbulk step, otherwise det can be in epsilon range
# Adjust config.yaml accordingly

#Run deconv
camObj <- CAM(T, K = length(rownames(P)), cores = 30)

#Extract results
res <- t(camObj@ASestResult[[1]]@Aest)

#Convert P to matrix
P <- matrix(unlist(P), nrow = nrow(P), ncol = ncol(P))

#Calculate RMSE error
rmse <- sqrt(mean((P - res)^2))

#Calculate euclidean distance (switch to any minkowski-type by adjusting p)
m_dist <- dist(rbind(as.vector(res), as.vector(unlist(P))), method = "minkowski", p = 2)

#Calculate Distance corr
distance_corr <- dcor(P, res)

#print and exit (update later)
print(paste0("RMSE: ", rmse))
print(paste0("Euclidean Distance: ", m_dist))
print(paste0("Distance Correlation: ", distance_corr))
