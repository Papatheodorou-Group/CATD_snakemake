## Pseudobulk generation script
##
## @zgr2788

library(Seurat)

#Read data
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1] #Name of T_ref
mode <- args[2] #Mode to be used
cellCount <- args[3] #How many cells to pool for each psuedobulk
nSamples <- args[4] #How many samples to generate
T_prep <- readRDS(filename)
filename <- sub("Input/Cell_splits", "Input/Psuedobulks", filename)

#set.seed(42)


#Switch modes
T_df <- data.frame(row.names = 1:nrow(T_prep@assays$RNA@counts))
P_df <- data.frame(row.names = levels(as.factor(T_prep@meta.data$cellType)))

switch(mode,
'1'={


  for (i in 1:nSamples)
  {
    #Get identifier for the sample
    identifierString <- paste0('sample_', i)

    #Classic random generation, pick cells (columns) randomly with replacement
    whichCols <- sample(seq_len(ncol(T_prep@assays$RNA@counts)), cellCount, replace = TRUE)

    #Get cell types corresponding to selected cells
    whichTypes <- T_prep@meta.data$cellType[whichCols]

    #Get picked columns and sum over them
    T <- rowSums(T_prep@assays$RNA@counts[,whichCols])

    #Get picked cell types and turn counts into proportions
    P <- table(whichTypes) / as.numeric(cellCount)

    #Append the bulk to dataframes
    T_df[identifierString] <- T
    P_df[identifierString] <- P


  }


},

'2'={
#Implement mode 2 here

}
)


#Debug - remove later
write.csv(T_df, "pbulks.csv")
write.csv(P_df, "props.csv")

#Save matrix to rds
saveRDS(T_df, file = sub("_gen.rds", "_pbulks.rds", filename))
saveRDS(P_df, file = sub("_gen.rds", "_props.rds", filename))
