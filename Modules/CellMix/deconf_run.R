## deconf deconv script
##
## @zgr2788

suppressMessages(library(CellMix))
####################################


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C2 <- args[3]
filename_O <- args[4]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C2 <- readRDS(filename_C2)


#Toss out the genes tossed out in T normalization from C as well
common <- intersect(rownames(C2), rownames(T))
C2 <- C2[common,]
T  <- T[common,]

#Preprocess
ML <- CellMix::MarkerList()
ML@.Data <- tapply(as.character(C2$geneName),as.character(C2$cellType),list)

#Get res and reorder the matrices for correspondence
res <- CellMix::ged(T, x=length(unique(as.character(C2$cellType))), method = "deconf", maxIter = 500, verbose= TRUE)
res <- CellMix::match.nmf(res, ML)@fit@H
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
