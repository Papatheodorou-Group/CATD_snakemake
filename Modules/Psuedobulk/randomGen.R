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



#Switch modes
switch(mode,
'1'={
  samplesDf <- data.frame(row.names = 1:nrow(T_prep@assays$RNA@counts))

  for (i in 1:nSamples)
  {
    #Get identifier for the sample
    identifierString <- paste0('sample_', i)

    #Classic random generation, pick cells (columns) randomly with replacement
    whichCols <- sample(seq_len(ncol(T_prep@assays$RNA@counts)), cellCount, replace = TRUE)

    #Get picked columns and sum over them
    T <- rowSums(T_prep@assays$RNA@counts[,whichCols])

    #Append the bulk to dataframe
    samplesDf[identifierString] <- T


  }

  #Debug - remove later
  write.csv(samplesDf, "samples.csv")


},

'2'={
#Implement mode 2 here

}
)
