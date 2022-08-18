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
filename_O <- args[4]


T <- readRDS(filename_T)
P <- readRDS(filename_P)


# Warning: Since fit model is linear, needs more samples in pbulk step, otherwise det can be in epsilon range
# Adjust config.yaml accordingly

#Run deconv
camObj <- CAM(T, K = length(rownames(P)), cores = as.numeric(args[5]))

#Extract results
res <- t(camObj@ASestResult[[1]]@Aest)

#Save and exit
saveRDS(res, file=filename_O)
