## MuSiC deconv script
##
## @zgr2788


#Load MuSiC
suppressMessages(library(devtools))
devtools::install_github('xuranw/MuSiC')
suppressMessages(library(MuSiC))
suppressMessages(library(Biobase))
suppressMessages(library(energy))
suppressMessages(library(dplyr))
suppressMessages(library(SeuratObject))



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
C0 <- ExpressionSet(C0, phenoData = as(phenData, "AnnotatedDataFrame"))


#Get results and reorder the matrices for correspondence
res <- t(music_prop(bulk.eset = T, sc.eset = C0, clusters = "cellType", samples = "sampleID", select.ct = levels(C0@phenoData$cellType), markers = NULL, normalize = FALSE, verbose = FALSE)$Est.prop.weighted)
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
