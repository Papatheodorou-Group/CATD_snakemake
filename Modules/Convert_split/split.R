## Split mini script in preparation for the Pseudobulk & C_gen modules
##
## @zgr2788

suppressMessages(library(Seurat))
#set.seed(42) #For reproducability purposes

partitionProp <- 0.5 #Can be edited later

filename <- commandArgs(trailingOnly = TRUE)
seuratObj <- readRDS(filename)
filename <- sub("Input", "Input/Cell_splits", filename)

splitSize <- floor(partitionProp * nrow(seuratObj@meta.data)) #Get the split size
idxList <- seq_len(nrow(seuratObj@meta.data)) #List of all indices

#Define low cell type counts
riskGroup <- names(table(seuratObj@meta.data$cellType)[table(seuratObj@meta.data$cellType) < 10])

#Assign split group to each cell randomly and split
seuratObj@meta.data$splitGroup <- as.integer(idxList %in% sample(idxList, size = splitSize))
seuratObj@meta.data$riskTag <- seuratObj@meta.data$cellType %in% riskGroup

splitRef <- seuratObj[,(seuratObj@meta.data$splitGroup == 1 | seuratObj@meta.data$riskTag == TRUE)]
splitGen <- seuratObj[,(seuratObj@meta.data$splitGroup == 0 | seuratObj@meta.data$riskTag == TRUE)]
#Check to see if all cells are sampled in one of the groups!


#Save the splits into outputs
saveRDS(splitGen, file = sub("_seurat.rds", "_gen.rds", filename)) #Sent to pseudobulk generation workflow
saveRDS(splitRef, file = sub("_seurat.rds", "_C0.rds", filename)) #Sent to ref/marker generation workflow
