#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C0 <- args[3]
filename_phenData <- args[4]
filename_O <- args[5]
filename_1 <- args[6]
T <- readRDS(filename_T)
P <- readRDS(filename_P)
C0 <- readRDS(filename_C0)
phenData <- readRDS(filename_phenData)


#Match genes in rows for both references
common <- intersect(rownames(C0), rownames(T))
T <- T[common,]
C0 <- C0[common,]
rm(common)

#Build sig
cellTypes <- phenData$cellType


# Assuming C0 is a dgCMatrix and phenData is a data.frame


C0 <- as.matrix(C0)


colnames(C0) <- cellTypes

C0<- rbind(colnames(C0), C0)
C0<- as.data.frame(C0)
C0$GeneSymbol <- rownames(C0)

C0 <- C0[, c("GeneSymbol", names(C0)[-ncol(C0)])]

C0[1,1]<-"GeneSymbol"



#Mixture T dataframe add a Gene column

Gene <- rownames(T)
T<-cbind(Gene,T)

#Make "Gene" the first column 
#T <- T[, c("Gene", names(T)[-ncol(T)])]

df <- T[, c(ncol(T), 1:(ncol(T)-1))]
#

write.table(C0, file = filename_O, sep="\t", col.names = FALSE, row.names = FALSE)

write.table(T, file = filename_1, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


