## Reference matrix type 2 (markergenes) generation script
##
## @zgr2788

suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(MAST))
suppressMessages(library(DESeq2))

plan('multisession', workers = 5) #Paralellism

#Read data
 args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
test_1 <- args[2]
test_2 <- args[3]
C_0 <- readRDS(filename)
filename <- sub("Input/Cell_splits", "Input/References", filename)

message(paste0('Running DGE with tests: ', test_1, ' & ', test_2))

#Port cell types to Idents
Idents(C_0) <- C_0@meta.data$cellType

#Find cell type specific markers
markers_t1 <- FindAllMarkers(C_0, min.pct = 0.5, logfc.threshold = log(2), test.use = test_1)
markers_t2 <- FindAllMarkers(C_0, min.pct = 0.5, logfc.threshold = log(2), test.use = test_2)

#Filter markers by corrected p-val
markers_t1 <- markers_t1[markers_t1$p_val_adj <= 0.05,]
markers_t2 <- markers_t2[markers_t2$p_val_adj <= 0.05,]

#Select those in both
markers <- markers_t1[rownames(markers_t1) %in% rownames(markers_t2),]

#Select cols
markers <- markers[,c("avg_log2FC", "cluster", "gene")]

#Change colnames
colnames(markers) <- c("log2FC", "cellType", "geneName")


#Debug - remove later
write.csv(markers, "C_2.csv")

#Write to rds
saveRDS(markers, file = sub("_C0", "_C2", filename ))
