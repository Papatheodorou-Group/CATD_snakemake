## TIMER deconv script
##
## @zgr2788


#Load TIMER
source("Modules/TIMER/TIMER.R")
suppressMessages(library(energy))
suppressMessages(library(Matrix))




#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C0 <- args[3]
filename_phenData <- args[4]
filename_O <- args[5]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C0 <- readRDS(filename_C0)
phenData <- readRDS(filename_phenData)


#Match genes in rows for both references
common <- intersect(rownames(C0), rownames(T))
T <- T[common,]
C0 <- C0[common,]
rm(common)

#Create ref_annotations
ref_anno <- phenData$cellID
names(ref_anno) <- phenData$cellType


#Get results and reorder the matrices for correspondence
res <-t(TIMER_deconv(mix = T, ref = C0, curated.cell.types = ref_anno, sig = rownames(T)))
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
