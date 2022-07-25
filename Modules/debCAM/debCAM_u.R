## debCAM unsupervised clustering script
##
## @zgr2788


suppressMessages(library(debCAM))
suppressMessages(library(BiocParallel))
register(MulticoreParam(4))
register(SnowParam(4))


#Read data
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
cellType_n <- as.numeric(args[2])
filename_P <- args[3]

T <- readRDS(filename_T)
P <- readRDS(filename_P)


# Warning: Since fit model is linear, needs at least n(cellTypes) samples in pbulk step, otherwise cannot get inverse
# Adjust config.yaml accordingly

#Run deconv
camObj <- CAM(T, K = cellType_n, thres.low = 0.3, thres.high = 0.95)

#Extract results
res <- t(camObj@ASestResult[[1]]@Aest)

#Fill empty rows if cell types were sampled
for (i in 1:(nrow(P)-nrow(res))) { res <- rbind(res, 0) }

#Calculate RMSE error
rmse <- sqrt(sum((P - res)^2) / length(res))

#Calculate euclidean distance (switch to any minkowski-type by adjusting p)
m_dist <- dist(rbind(as.vector(res), as.vector(unlist(P))), method = "minkowski", p = 2)

#print and exit (update later)
print(paste0("RMSE: ", rmse))
print(paste0("Euclidean Distance: ", m_dist))
