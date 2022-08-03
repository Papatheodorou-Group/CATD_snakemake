## debCAM semi-supervised deconv script with cell type profile matrix (S in debCAM documentation)
##
## @zgr2788


suppressMessages(library(debCAM))
suppressMessages(library(energy))


#Read data
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1 <- args[3]
filename_O <- args[4]

T <- readRDS(filename_T)
C1 <- readRDS(filename_C1)
P <- readRDS(filename_P)

#Match cols in C1 to rownames in P for true reference
C1 <- C1[,order(match(colnames(C1), rownames(P)))]

#Get foldchanges for possible marker genes
pMGstat <- MGstatistic(C1, colnames(C1))

#Get those with a fold change over 10, might need to revisit
pMGlist.FC <- lapply(colnames(C1), function(x) rownames(pMGstat)[pMGstat$idx == x & pMGstat$OVE.FC > 10])

#Run with marker list generated from S
res <- redoASest(T, pMGlist.FC, maxIter = 10) #Adjust maxIter later
res <- t(res$Aest)

#Save and exit
saveRDS(res, file=filename_O)
