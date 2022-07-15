## Split mini script in preparation for the Pseudobulk & C_gen modules
##
## @zgr2788

library(Seurat)
set.seed(42) #For reproducability purposes

partitionProp <- 0.5 #Can be edited later

filename <- commandArgs(trailingOnly = TRUE)
seuratObj <- readRDS(filename)
filename <- sub("Scratch/Convert_split", "Output/Cell_splits", filename)

splitSize <- floor(partitionProp * nrow(seuratObj@meta.data)) #Get the split size
idxList <- seq_len(nrow(seuratObj@meta.data)) #List of all indices

#Assign split group to each cell randomly and split
seuratObj@meta.data$splitGroup <- as.integer(idxList %in% sample(idxList, size = splitSize))
seuratSplit <- SplitObject(seuratObj, split.by = 'splitGroup')

#Save the splits into outputs
saveRDS(seuratSplit[[1]], file = sub(".rds", "_gen.rds", filename)) #Sent to pseudobulk generation workflow
saveRDS(seuratSplit[[2]], file = sub(".rds", "_C0.rds", filename)) #Sent to ref/marker generation workflow
