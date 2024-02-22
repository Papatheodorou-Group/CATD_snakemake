## BayesPrism deconv script
##
## @zgr2788
## @annavpo


#Load omnideconv
#Make sure BayesPrism dependency is installed. 
suppressMessages(library(omnideconv))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C0 <- args[3]
filename_phenData <- args[4]
cores <- as.numeric(args[5])
filename_O <- args[6]

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
res <- omnideconv::deconvolute(T, signature = NULL, cell_type_annotations=cellTypes, single_cell_object = C0, method = "bayesprism",n_cores=cores)

res[res < 10^-5] <- 0 #Convergence error tolerance = 10^-5
res <- t(res/rowSums(res))

if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
