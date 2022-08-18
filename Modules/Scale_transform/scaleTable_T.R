## Matrix scaling script for T  matrices
##
## @zgr2788

suppressMessages(library(Matrix))
suppressMessages(library(future))


#Read data
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
method <- args[2]
plan('multisession', workers = as.numeric(args[3])) #Paralellism


T <- readRDS(filename_T)
filename_T <- sub("Input/Psuedobulks", "Input/Normalized_tables", filename_T)


#Preprocess matrices to avoid errorshaving multiple snakefiles in one directory

message(paste0("Scaling T with method: ", method))

#Only get rows where there is at least one cell with a count
T <- T[rowSums(T) != 0,]

#Same for cols
T <- T[,colSums(T) != 0]

#Only keep rows with different counts after transformations
T <- T[!apply(T, 1, function(x) var(x) == 0),]


#Switch modes
switch(method,


  "none" = {
    message('Warning: T not scaled')
  },


  "column" = {
    T <- apply(T, 2, function(x) x/sum(x))
  },


  "row" = {
    T <- t(apply(T, 1, function(x) x/sum(x)))
  },


  "mean" = {
    T <- apply(T, 2, function(x) x - mean(x))
  },


  "column_z-score" = {
    T <- scale(T, center = TRUE, scale = TRUE)
  },


  "global_z-score" = {
    T <- (as.matrix(T) - mean(as.matrix(T))) / sd(as.matrix(T))
  },


  "column_min-max" = {
    T <- apply(T, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  },


  "global_min-max" = {
    T <- (T - min(T))/(max(T) - min(T))
  },


  "LogNormalize" = {
    suppressMessages(library(Seurat))

    T <- as.matrix(expm1(Seurat::LogNormalize(T, verbose = FALSE)))
  },


  "QN" = {
    suppressMessages(library(preprocessCore))

    T_rownames <- rownames(T); T_colnames <- colnames(T)
    T <- preprocessCore::normalize.quantiles(as.matrix(T))
    rownames(T) <- T_rownames; colnames(T) <- T_colnames
  },


  "TMM" = {
    suppressMessages(library(edgeR))

    T <- edgeR::DGEList(counts = T, group = colnames(T))
    T <- edgeR::calcNormFactors(T, method = "TMM")
    T <- edgeR::cpm(T)
  },


  "UQ" = {
    suppressMessages(library(edgeR))

    T <- edgeR::DGEList(counts = T, group = colnames(T))
    T <- edgeR::calcNormFactors(T, method = "upperquartile")
    T <- edgeR::cpm(T)
  },


  "median_ratios" = {
    suppressMessages(library(DESeq2))

    T <- round(T); metadata <- colnames(T);  #DESeq2 requires integers as counts
    dds <- DESeq2::DESeqDataSetFromMatrix(T, colData = metadata, design = ~samples)
    dds <- DESeq2::estimateSizeFactors(dds, type = "ratio")
    T <- DESeq2::counts(dds, normalized = TRUE)
  },


#  "TPM" = {
#    suppressMessages(library(SingleR))
#
#    #Cannot find data(human_lengths) ? -> check later
#  }


)

#Only get rows where there is at least one cell with a count
T <- T[rowSums(T) != 0,]

#Same for cols
T <- T[,colSums(T) != 0]

#Only keep rows with different counts after transformations
T <- T[!apply(T, 1, function(x) var(x) == 0),]

#Save transformed tables
saveRDS(T, sub(".rds", "_scaled.rds",filename_T))

#Debug
#write.csv(T, "T_scaled.csv")
