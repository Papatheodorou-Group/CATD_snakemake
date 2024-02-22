## proportionsInAdmixture deconv script
##
## @zgr2788


#Load ADAPTS
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("preprocessCore", configure.args = c(preprocessCore = "--disable-threading"), force= TRUE, update=TRUE, type = "source")
if (!require("ADAPTS", quietly = TRUE))
    install.packages("ADAPTS", repos='http://cran.us.r-project.org')
if (!require("WGCNA", quietly = TRUE))
    install.packages("WGCNA", repos='http://cran.us.r-project.org')


suppressMessages(library(ADAPTS))
suppressMessages(library(WGCNA))
#suppressMessages(library(energy))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
filename_P <- args[2]
filename_C1<- args[3]
filename_C2 <- args[4]
markers <- as.logical(args[5])
filename_O <- args[6]

T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)
C2 <- readRDS(filename_C2)


if (markers) { C1 <- C1[C2$geneName,]}


#Toss out the genes tossed out in T normalization from C as well
common <- intersect(rownames(C1), rownames(T))
C1 <- C1[common,]
T  <- T[common,]


#Get res and reorder the matrices for correspondence
res <- ADAPTS::estCellPercent(refExpr = C1, geneExpr = T, method="proportionsInAdmixture")
res[is.na(res)] <- 0
res = apply(res,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
res = apply(res,2,function(x) x/sum(x)) #explicit STO constraint
res <- res[-nrow(res),]
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)
