## CDSeq deconv script
##
## @zgr2788


#Load ADAPTS
if (!require("ADAPTS", quietly = TRUE))
    install.packages("ADAPTS", repos='http://cran.us.r-project.org')
if (!require("WGCNA", quietly = TRUE))
    install.packages("WGCNA", repos='http://cran.us.r-project.org')


suppressMessages(library(ADAPTS))
suppressMessages(library(WGCNA))
suppressMessages(library(energy))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1 <- args[3]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)


#Toss out the genes tossed out in T normalization from C as well
C1 <- C1[rownames(C1) %in% rownames(T),]

#Get res and reorder the matrices for correspondence
res <- ADAPTS::estCellPercent(refExpr = C1, geneExpr = T, method="proportionsInAdmixture")
res[is.na(res)] <- 0
res = apply(res,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
res = apply(res,2,function(x) x/sum(x)) #explicit STO constraint
res <- res[-nrow(res),]
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
