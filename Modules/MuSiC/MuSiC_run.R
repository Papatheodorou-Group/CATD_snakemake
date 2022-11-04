## MuSiC deconv script
##
## @zgr2788
library('devtools')

#install SingleCellExperiment
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment",force=TRUE)

#install TOAST
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TOAST",force=TRUE)

#install Biobase
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biobase",force=TRUE)

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')

install.packages('energy',repos='http://cran.us.r-project.org')
suppressMessages(library(MuSiC))
suppressMessages(library(Biobase))
suppressMessages(library(energy))
suppressMessages(library(dplyr))
install.packages('SeuratObject',repos='http://cran.us.r-project.org')
suppressMessages(library(SeuratObject))
suppressMessages(library(SingleCellExperiment))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C0 <- args[3]
filename_phenData <- args[4]
force_raw <- as.logical(args[5])
filename_O <- args[6]

if (!force_raw)
  {
    message("WARNING: MuSiC requires RAW counts to run correctly.\nIf you did not intend to run MuSiC with scaling/transformations applied, set\n\n'force_raw = TRUE'\n\nin config.yaml.")
  } else {
    message("MuSiC is running with RAW counts...")
    filename_T <- filename_T %>% sub("Normalized_tables", "Psuedobulks" , .) %>% sub("_scaled_transformed", "", .) %>% sub("_transformed_scaled", "", .)
    filename_C0 <- filename_C0 %>% sub("Normalized_tables", "Cell_splits", .) %>% sub("_scaled_transformed", "", .) %>% sub("_transformed_scaled", "", .)
  }

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C0 <- readRDS(filename_C0)
phenData <- readRDS(filename_phenData)
if (force_raw) { C0 <- as.matrix(C0@assays$RNA@counts) ; T <- as.matrix(T) }

#Match genes in rows for both references
common <- intersect(rownames(C0), rownames(T))
T <- T[common,]
C0 <- C0[common,]
rm(common)

#Convert to eset
T <- ExpressionSet(T)
print(str(T))
C0 <- SingleCellExperiment(list(counts=as.matrix(C0)),colData=DataFrame(cellType=phenData$cellType,sampleID=phenData$sampleID),
    metadata=list(phenData))
print(str(C0))
T <- exprs(T)
#Get results and reorder the matrices for correspondence
res <- music_prop(T, C0, clusters = "cellType", samples = "sampleID", markers = NULL, normalize = FALSE, verbose = TRUE)$Est.prop.weighted
res <- t(res)
print(str(res))
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
