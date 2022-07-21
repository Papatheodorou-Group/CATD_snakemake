## Reference matrix type 1 (rows = genenames, cols = celltypes) generation script
##
## @zgr2788

suppressMessages(library(Seurat))
suppressMessages(library(sparseMatrixStats))

#Read data
filename <- commandArgs(trailingOnly = TRUE)
C_0 <- readRDS(filename)
filename <- sub("Input/Cell_splits", "Input/References", filename)

#Split cells by type
cellSplits <- SplitObject(C_0, split.by = "cellType")

#Initialize references
C_1 <- data.frame(row.names = rownames(C_0@assays$RNA@counts))
refVar <- data.frame(row.names = rownames(C_0@assays$RNA@counts))

#Group all celltype rowsums
for (i in 1:length(levels(C_0@meta.data$cellType)))
{
    C_1[names(cellSplits[i])] <- rowMeans(cellSplits[[i]]@assays$RNA@counts)
    refVar[names(cellSplits[i])] <- sparseMatrixStats::rowSds(cellSplits[[i]]@assays$RNA@counts)
}

#Write sc-counts matrix as separate object to save memory later on
C_counts <- C_0@assays$RNA@counts

#Debug
write.csv(C_1, "C_1.csv")
write.csv(refVar, "C_refVar.csv")

#Write to RDS
saveRDS(C_1, file = sub("_C0", "_C1", filename ))
saveRDS(refVar, file = sub("_C0", "_refVar", filename ))
saveRDS(C_counts, file = filename)
