## Matrix transformation script for T & C matrices
##
## @zgr2788

suppressMessages(library(DESeq2))
suppressMessages(library(future))



#Read data
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_C <- args[2]
method <- args[3]
plan('multisession', workers = as.numeric(args[4])) #Paralellism


T <- readRDS(filename_T)
C <- readRDS(filename_C)

#Debug for local, remove later
#C <- C[,1:500]

#Correction for filenames if needed
filename_T <- sub("Input/Psuedobulks", "Input/Normalized_tables", filename_T)
filename_C <- sub("Input/References", "Input/Normalized_tables", filename_C)


message(paste0("Transforming tables with method: ", method))

switch(method,

  "none" = {
    message('Warning: No transformation applied')
  },

  "log" = {
    T <- log1p(T)
    C <- log1p(C)
  },

  "sqrt" = {
    T <- sqrt(T)
    C <- sqrt(C)
  },

  "vst" = {
    #Not recommended for low memory, uses dense matrix
    #vst requires int, convert first
    T <- round(T)
    C <- round(C)

    T <- DESeq2::varianceStabilizingTransformation(T)
    C <- DESeq2::varianceStabilizingTransformation(C)
  }

)

#Only get rows where there is at least one cell with a count
C <- C[rowSums(C) != 0,]

#Same for cols
C <- C[,colSums(C) != 0]

#Only keep rows with different counts after transformations
C <- C[!apply(C, 1, function(x) var(x) == 0),]

#Only get rows where there is at least one cell with a count
T <- T[rowSums(T) != 0,]

#Same for cols
T <- T[,colSums(T) != 0]

#Only keep rows with different counts after transformations
T <- T[!apply(T, 1, function(x) var(x) == 0),]

#Save transformed tables
saveRDS(T, sub(".rds", "_transformed.rds", filename_T))
saveRDS(C, sub(".rds", "_transformed.rds", filename_C))

#Debug
#write.csv(T, "T_transformed.csv")
