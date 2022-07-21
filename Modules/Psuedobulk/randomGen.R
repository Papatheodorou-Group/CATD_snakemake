## Pseudobulk generation script
##
## @zgr2788

suppressMessages(library(Seurat))

#Read data
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1] #Name of T_ref
mode <- args[2] #Mode to be used
cellCount <- as.numeric(args[3]) #How many cells to pool for each psuedobulk
nSamples <- as.numeric(args[4]) #How many samples to generate
tryCatch(
  expr={
    propVar <- as.numeric(args[5])
    sampleCT <- as.numeric(args[6])
  }
  ,
  warning= function(warn){
    message("Warning: less than expected args for proportional sampling, running in random sampling...")
    mode <<- '1' #Altering global variable mode here
  }
)

T_prep <- readRDS(filename)
filename <- sub("Input/Cell_splits", "Input/Psuedobulks", filename)

#set.seed(42)


#Switch modes
T_df <- data.frame(row.names = rownames(T_prep@assays$RNA@counts))
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

    for (i in 1:nSamples)
    {
      #Get identifier for the sample
      identifierString <- paste0('sample_', i)



      #Sampling cell types
      if (sampleCT)
      {
        includedTypeCount <- as.numeric(sample(3:nlevels(as.factor(T_prep@meta.data$cellType)), 1)) #Minimum 3 cell types needed
        toKeep <- sample(levels(T_prep@meta.data$cellType), includedTypeCount, replace = FALSE) #Sample the cell types
        T_prep@meta.data$toKeep <- as.numeric(T_prep@meta.data$cellType %in% toKeep) #Assign kept/unkept to cells
        keepObj <- SplitObject(T_prep, split.by = 'toKeep') #Keep the selected cell types
        if(all(keepObj[[1]]@meta.data$toKeep == 1)) { keepObj <- keepObj[[1]] } else { keepObj <- keepObj[[2]] } #Need to differentiate between kept cell types
        keepObj@meta.data$cellType <- droplevels(keepObj@meta.data$cellType) #Drop unused cell types
      }

      else { keepObj <- T_prep }



      #Applying distributions
      if (propVar <= 0)
      {
        #Sample with min max bound 1-99
        P <- runif(nlevels(keepObj@meta.data$cellType), 1, 99)
      }

      else
      {
        #Fit the distribution with given variance
        P <- round(runif(nlevels(keepObj@meta.data$cellType), 100 - propVar, 100 + propVar))/100
      }

      #Sum to 1
      P <- round(P/sum(P), digits = log10(cellCount))

      #Get how many cells will be sampled by multiplying with prop matrix
      correspondingCounts <- cellCount * P

      #Sample to create T
      splitObj <- SplitObject(keepObj, split.by = 'cellType') #Split by cell type
      T <- numeric(length = nrow(keepObj@assays$RNA@counts)) #Initialize T

      #Add samples into T
      for (j in 1:length(correspondingCounts))
      {
        #Get sample per cell type
        whichCols <- sample(seq_len(ncol(splitObj[[j]]@assays$RNA@counts)), correspondingCounts[[j]], replace = TRUE)
        T <- T + rowSums(splitObj[[j]]@assays$RNA@counts[,whichCols])
        P_df[names(splitObj[j]), identifierString] <- P[[j]]
      }

      #Append the bulk to dataframes
      T_df[identifierString] <- T
      P_df[is.na(P_df)] <- 0
    }

  }
)


#Debug - remove later
write.csv(T_df, "pbulks.csv")
write.csv(P_df, "props.csv")

#Save matrix to rds
saveRDS(T_df, file = sub("_gen.rds", "_pbulks.rds", filename))
saveRDS(P_df, file = sub("_gen.rds", "_props.rds", filename))
