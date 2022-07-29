## DeconRNASeq deconv script
##
## @zgr2788


#Load DeconRNASeq
suppressMessages(library(BiocManager))
BiocManager::install('DeconRNASeq')
suppressMessages(library(energy))
suppressMessages(library(DeconRNASeq))



#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1 <- args[3]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)

#Explicit typecast to dataframe
T <- data.frame(T)
C1 <- data.frame(C1)

#Run DeconRNASeq and reorder rows for correspondence
res <- DeconRNASeq(datasets = T, signatures = C1, proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = TRUE, fig = FALSE)
res <- t(res$out.all)
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