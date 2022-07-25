## Reference matrix type 1 (rows = genenames, cols = celltypes) generation script
##
## @zgr2788

suppressMessages(library(Seurat))
suppressMessages(library(sparseMatrixStats))

#Read data
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
transform <- args[2]
C_0 <- readRDS(filename)
filename <- sub("Input/Cell_splits", "Input/References", filename)

#Write sc-counts matrix as separate object to save memory later on
C_counts <- C_0@assays$RNA@counts
C_metadata <- C_0@meta.data

#Initialize references
C_1 <- data.frame(row.names = rownames(C_0@assays$RNA@counts))
refVar <- data.frame(row.names = rownames(C_0@assays$RNA@counts))

#Transform C matrix
switch(transform,

  "none" = {
    message('Warning: No transformation applied')
  },

  "log" = {
    C_0@assays$RNA@counts <- log1p(C_0@assays$RNA@counts)
  },

  "sqrt" = {
    C_0@assays$RNA@counts <- sqrt(C_0@assays$RNA@counts)
  },

  "vst" = { #DO NOT USE
    #Not recommended for low memory, uses dense matrix
    #vst requires int, convert first
    C_0@assays$RNA@counts <- round(C_0@assays$RNA@counts)

    C_0@assays$RNA@counts <- DESeq2::varianceStabilizingTransformation(C_0@assays$RNA@counts)
  }

)

#Split cells by type
cellSplits <- SplitObject(C_0, split.by = "cellType")

#Group all celltype rowsums
for (i in 1:length(levels(C_0@meta.data$cellType)))
{
    C_1[names(cellSplits[i])] <- rowMeans(cellSplits[[i]]@assays$RNA@counts)
    refVar[names(cellSplits[i])] <- sparseMatrixStats::rowSds(cellSplits[[i]]@assays$RNA@counts)
}

#Debug
#write.csv(C_1, "C_1.csv")
#write.csv(refVar, "C_refVar.csv")

#Write to RDS
saveRDS(C_1, file = sub("_C0", "_C1", filename ))
saveRDS(refVar, file = sub("_C0", "_refVar", filename ))
saveRDS(C_counts, file = filename)
saveRDS(C_metadata, file = sub("_C0", "_phenData", filename))
