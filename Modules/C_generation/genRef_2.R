## Reference matrix type 2 (markergenes) generation script
##
## @zgr2788

library(Seurat)

#Read data
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
topMarkers <- as.numeric(args[2]) #Useless as of now
C_0 <- readRDS(filename)
filename <- sub("Input/Cell_splits", "Input/References", filename)

#Port cell types to Idents
Idents(C_0) <- C_0@meta.data$cellType

#Do DGE analysis on marker groups

#Find cell type specific markers
markers <- FindAllMarkers(C_0, min.pct = 0.5, logfc.threshold = log(2), test.use = "wilcox")

#Select cols
markers <- markers[,c("avg_log2FC", "cluster", "gene")]

#Change colnames
colnames(markers) <- c("log2FC", "cellType", "geneName")


#Debug - remove later
write.csv(markers, "C_2.csv")

#Write to rds
saveRDS(markers, file = sub("_C0.rds", "_C2.rds", filename ))
