## Matrix scaling script for C matrices
##
## @zgr2788

suppressMessages(library(Matrix))

#Read data
args <- commandArgs(trailingOnly = TRUE)
filename_C <- args[1]
method <- args[2]

C <- readRDS(filename_C)
filename_C<- sub("Input/References", "Input/Normalized_tables", filename_C)


#debug on local to reduce dimensions, delete later
C <- C[,1:500]


#Preprocess matrices to avoid errors

#Only get rows where there is at least one cell with a count
C <- C[rowSums(C) != 0,]

#Same for cols
C <- C[,colSums(C) != 0]

#Only keep rows with different counts after transformations
C <- C[!apply(C, 1, function(x) var(x) == 0),]


#Switch modes
switch(method,


  "none" = {
    #Do nothing
  },


  "column" = {
    C <- apply(C, 2, function(x) x/sum(x))
  },


  "row" = {
    C <- t(apply(C, 1, function(x) x/sum(x)))
  },


  "mean" = {
    C <- apply(C, 2, function(x) x - mean(x))
  },


  "column_z-score" = {
    C <- scale(C, center = TRUE, scale = TRUE)
  },


  "global_z-score" = {
    C <- (as.matrix(C) - mean(as.matrix(C))) / sd(as.matrix(C))
  },


  "column_min-max" = {
    C <- apply(C, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  },


  "global_min-max" = {
    C <- (C - min(C))/(max(C) - min(C))
  },


  "LogNormalize" = {
    suppressMessages(library(Seurat))

    C <- as.matrix(expm1(Seurat::LogNormalize(C, verbose = FALSE)))
  },


  "QN" = {
    suppressMessages(library(preprocessCore))

    C_rownames <- rownames(C); C_colnames <- colnames(C)
    C <- preprocessCore::normalize.quantiles(as.matrix(C))
    rownames(C) <- C_rownames; colnames(C) <- C_colnames
  },


  "CMM" = { #CHECK FOR C
    suppressMessages(library(edgeR))

    C <- edgeR::DGEList(counts = C, group = colnames(C))
    C <- edgeR::calcNormFactors(C, method = "TMM")
    C <- edgeR::cpm(C)
  },


  "UQ" = { #CHECK FOR C
    suppressMessages(library(edgeR))

    C <- edgeR::DGEList(counts = C, group = colnames(C))
    C <- edgeR::calcNormFactors(C, method = "upperquartile")
    C <- edgeR::cpm(C)
  },


  "median_ratios" = { #CHECK FOR C
    suppressMessages(library(DESeq2))

    C <- round(C); metadata <- colnames(C);  #DESeq2 requires integers as counts
    dds <- DESeq2::DESeqDataSetFromMatrix(C, colData = metadata, design = ~samples)
    dds <- DESeq2::estimateSizeFactors(dds, type = "ratio")
    C <- DESeq2::counts(dds, normalized = TRUE)
  },


  "TPM" = {
    suppressMessages(library(SingleR))

    #Cannot find data(human_lengths) ? -> check later
  }


)

#Save transformed tables
saveRDS(C, sub(".rds", "_scaled.rds",filename_C))
