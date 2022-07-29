## CIBERSORT deconv script
##
## @zgr2788


#Load CIBERSORT and prereqs
source('Modules/CIBERSORT/CIBERSORT.R')
suppressMessages(library(e1071))
suppressMessages(library(preprocessCore))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(energy))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1 <- args[3]


T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)



#Explicit typecasting for T to run with tibble
T <- data.frame(T)

#Run deconv, params may be exported to config later
res <- CIBERSORT(sig_matrix = C1, mixture_file = T, QN = FALSE, absolute = FALSE, perm = 0)
res <- t(res[,1:(ncol(res)-3)])
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
