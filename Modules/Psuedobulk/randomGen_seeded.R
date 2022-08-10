## Pseudobulk generation script
##
## @zgr2788

if (!require("FamilyRank", quietly = TRUE)){
    install.packages("FamilyRank", repos='http://cran.us.r-project.org')
}

suppressMessages(library(Seurat))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(dplyr))
suppressMessages(library(FamilyRank))






#Function used for parallel combinations
combFunc <- function(...) {
    mapply('bind_cols', ..., SIMPLIFY=FALSE)
}

#Read data
args <- commandArgs(trailingOnly = TRUE)

cores <- as.numeric(args[7])
registerDoParallel(cores)

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

set.seed(42)

T_prep <- readRDS(filename)
filename <- sub("Input/Cell_splits", "Input/Psuedobulks", filename)

#set.seed(42)


#Switch modes

switch(mode,
'1'={


    T_P <- foreach (i=1:nSamples, .combine = combFunc) %dopar%
    {

      #Classic random generation, pick cells (columns) randomly with replacement
      whichCols <- sample(seq_len(ncol(T_prep@assays$RNA@counts)), cellCount, replace = TRUE)

      #Get cell types corresponding to selected cells
      whichTypes <- T_prep@meta.data$cellType[whichCols]

      #Get picked columns and sum over them
      T <- data.frame(rowSums(T_prep@assays$RNA@counts[,whichCols]))
      P <- data.frame(table(whichTypes) / as.numeric(cellCount), row.names = 'whichTypes')
      list(T, P)
    }
    colnames(T_P[[1]]) <- paste0("sample_", 1:ncol(T_P[[1]]))
    colnames(T_P[[2]]) <- paste0("sample_", 1:ncol(T_P[[2]]))


},

'2'={

    T_P <- foreach (i=1:nSamples, .combine = combFunc) %dopar%
    {



      #Applying distributions
      if (sampleCT && (propVar > 0))
      {
        P <- abs(rbinorm(nlevels(T_prep@meta.data$cellType), 0, 1, 0.0001, propVar, 0.5))
        P <- as.numeric(lapply(P, function(x) ifelse(x == 0, 1e-9, x)))
      }

      else if (propVar <= 0)
      {
        #Sample with min max bound 1-99
        P <- runif(nlevels(T_prep@meta.data$cellType), 1, 99)
      }

      else
      {
        #Fit the distribution with given variance
        P <- round(runif(nlevels(T_prep@meta.data$cellType), 100, 100 + (propVar*2)))/100
      }

      #Sum to 1
      P <- P/sum(P)

      #Get how many cells will be sampled by multiplying with prop matrix
      correspondingCounts <- round(cellCount * P)

      #Sample to create T
      splitObj <- SplitObject(T_prep, split.by = 'cellType') #Split by cell type
      T <- numeric(length = nrow(T_prep@assays$RNA@counts)) #Initialize T
      P <- data.frame(P, row.names = names(splitObj))

      #Add samples into T
      for (j in 1:length(correspondingCounts))
      {
        #Get sample per cell type
        whichCols <- sample(seq_len(ncol(splitObj[[j]]@assays$RNA@counts)), correspondingCounts[[j]], replace = TRUE)
        tryCatch(
          expr={
            T <- T + rowSums(splitObj[[j]]@assays$RNA@counts[,whichCols])
          }
          ,
          error= function(e){
          }
        )
      }

      #Append the bulk to dataframes
      T <- data.frame(T)

      list(T, P)
    }
    colnames(T_P[[1]]) <- paste0("sample_", 1:ncol(T_P[[1]]))
    colnames(T_P[[2]]) <- paste0("sample_", 1:ncol(T_P[[2]]))

  }
)


#Debug - remove later
#write.csv(T_df, "pbulks.csv")
#write.csv(P_df, "props.csv")

#Save matrix to rds
saveRDS(T_P[[1]], file = sub("_gen_seeded.rds", "_pbulks_seeded.rds", filename))
saveRDS(T_P[[2]], file = sub("_gen_seeded.rds", "_props_seeded.rds", filename))
