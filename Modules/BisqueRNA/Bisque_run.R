## Bisque deconv script
##
## @zgr2788


#Load Bisque
suppressMessages(library(BisqueRNA))
suppressMessages(library(Biobase))
suppressMessages(library(energy))



#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C0 <- args[3]
filename_phenData <- args[4]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C0 <- readRDS(filename_C0)
phenData <- readRDS(filename_phenData)


#Match genes in rows for both references
common <- intersect(rownames(C0), rownames(T))
T <- T[common,]
C0 <- C0[common,]
rm(common)

#Convert to eset
T <- ExpressionSet(T)
C0 <- ExpressionSet(C0, phenoData = as(phenData, "AnnotatedDataFrame"))
C0$SubjectName <- C0$sampleID #BisqueRNA uses hard-coded colname


#Get results and reorder the matrices for correspondence
res <- BisqueRNA::ReferenceBasedDecomposition(T, C0, markers=NULL, use.overlap=FALSE)$bulk.props
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
