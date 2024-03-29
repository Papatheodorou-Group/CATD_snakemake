## CDSeq deconv script
##
## @zgr2788


#Load CDseq
suppressMessages(library(devtools))
devtools::install_github("kkang7/CDSeq_R_Package")
suppressMessages(library(CDSeq))
suppressMessages(library(energy))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1 <- args[3]
filename_C2 <- args[4]
cores <- as.numeric(args[5])
markers <- as.logical(args[6])
filename_O <- args[7]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)
C2 <- readRDS(filename_C2)

if (markers) { C1 <- C1[C2$geneName,]}

#Toss out the genes tossed out in T normalization from C as well
common <- intersect(rownames(C1), rownames(T))
C1 <- C1[common,]
T  <- T[common,]

#Get results and reorder the matrices for correspondence
res <- CDSeq(bulk_data = T, reference_gep = C1, cell_type_number = ncol(C1), mcmc_iterations = 1000, cpu_number = cores, block_number = 6, gene_subset_size=15)
res <- res$estProp
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
