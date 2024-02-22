## DeconRNASeq deconv script
##
## @zgr2788


#Load DeconRNASeq
suppressMessages(library(BiocManager))
BiocManager::install('DeconRNASeq')
suppressMessages(library(energy))
suppressMessages(library(DeconRNASeq))



#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1 <- args[3]
filename_C2 <- args[4]
markers <- as.logical(args[5])
filename_O <- args[6]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)
C2 <- readRDS(filename_C2)


if (markers) { C1 <- C1[C2$geneName,]}
#Explicit typecast to dataframe
T <- data.frame(T)
C1 <- data.frame(C1)

#Run DeconRNASeq and reorder rows for correspondence
res <- DeconRNASeq(datasets = T, signatures = C1, proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = TRUE, fig = FALSE)
res <- t(res$out.all)
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
