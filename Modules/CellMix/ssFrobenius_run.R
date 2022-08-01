## CDSeq deconv script
##
## @zgr2788


#Load CellMix
if (!require("CellMix", quietly = TRUE)){
    install.packages("Modules/CellMix/BiocInstaller_1.8.3.tar.gz", repos=NULL, type = "source")
    install.packages("Modules/CellMix/CellMix_1.6.2_R_x86_64-conda-linux-gnu.tar.gz", repos=NULL, type = "source")
}

suppressMessages(library(CellMix))
suppressMessages(library(energy))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C2 <- args[3]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C2 <- readRDS(filename_C2)


#Toss out the genes tossed out in T normalization from C as well
C2 <- C2[rownames(C2) %in% rownames(T),]

#Preprocess
ML <- CellMix::MarkerList()
ML@.Data <- tapply(as.character(C2$geneName),as.character(C2$cellType),list)

#Get res and reorder the matrices for correspondence
res <- CellMix::ged(as.matrix(T), ML, method = "ssFrobenius", sscale = TRUE, maxIter = 500, log = FALSE)@fit@H 
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
