## Reference matrix type 1 (rows = genenames, cols = celltypes) generation script
##
## @zgr2788

#Read data
library(Seurat)
filename <- commandArgs(trailingOnly = TRUE)
C_0 <- readRDS(filename)
filename <- sub("Input/Cell_splits", "Input/References", filename)

#Split cells by type
cellSplits <- SplitObject(C_0, split.by = "cellType")

#Initialize reference
C_1 <- data.frame(row.names = rownames(C_0@assays$RNA@counts))

#Group all celltype rowsums
for (i in 1:length(levels(C_0@meta.data$cellType)))
{
    C_1[names(cellSplits[i])] <- rowMeans(cellSplits[[i]]@assays$RNA@counts)
}

#Debug
write.csv(C_1, "C_1.csv")

#Write to RDS
saveRDS(C_1, file = sub("_C0.rds", "_C1.rds", filename ))
