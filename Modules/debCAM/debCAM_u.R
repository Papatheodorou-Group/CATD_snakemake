## debCAM unsupervised clustering script
##
## @zgr2788


suppressMessages(library(debCAM))
suppressMessages(library(BiocParallel))


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
camObj <- CAM(T, K = cellType_n)

#Extract results
res <- t(camObj@ASestResult[[1]]@Aest)

#Fill empty rows if cell types were sampled
for (i in 1:(nrow(P)-nrow(res))) { res <- rbind(res, 0) }

#Convert P to matrix
P <- matrix(unlist(P), nrow = nrow(P), ncol = ncol(P))

#Calculate RMSE error
rmse <- sqrt(mean((P - res)^2))

#Calculate euclidean distance (switch to any minkowski-type by adjusting p)
m_dist <- dist(rbind(as.vector(res), as.vector(unlist(P))), method = "minkowski", p = 2)

#print and exit (update later)
print(paste0("RMSE: ", rmse))
print(paste0("Euclidean Distance: ", m_dist))
