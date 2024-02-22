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
filename_C2 <- args[5]
markers <- as.logical(args[6])
filename_O <- args[7]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C0 <- readRDS(filename_C0)
C2 <- readRDS(filename_C2)
phenData <- readRDS(filename_phenData)

#dense matrix
C0<-as.matrix(C0)
  



#Match genes in rows for both references
common <- intersect(rownames(C0), rownames(T))
T <- T[common,]
C0 <- C0[common,]
rm(common)

#Create ref_annotations
ref_anno <- phenData$cellID
names(ref_anno) <- phenData$cellType

if (markers) {
  sig <- intersect(C2$geneName, rownames(T))
} else {
  sig <- rownames(T)
}



#Get results and reorder the matrices for correspondence
res <-TIMER_deconv(mix = T, ref = C0, curated.cell.types = ref_anno, sig = sig)
res<- t(res)
res<-as.data.frame(res)
res = apply(res,2,function(x) ifelse(x < 0, 0, x)) #Explicit non-negativity constraint
res = apply(res,2,function(x) x/sum(x)) #Explicit STO constraint
res[is.na(res)] <- 0     

if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)


