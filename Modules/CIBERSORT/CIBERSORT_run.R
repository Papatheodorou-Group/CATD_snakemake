## CIBERSORT deconv script
##
## @zgr2788


#Load CIBERSORT and prereqs
source('Modules/CIBERSORT/CIBERSORT.R')
suppressMessages(library(e1071))
suppressMessages(library(preprocessCore))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(energy))


#Get args and load files
args <- commandArgs(trailingOnly = TRUE)
filename_T <- args[1]
print(filename_T)
filename_P <- args[2]
print(filename_P)
filename_C1 <- args[3]
print(filename_C1)
filename_C2 <- args[4]
print(filename_C2)
markers <- as.logical(args[5])
print(markers)
filename_O <- args[6]


T <- readRDS(filename_T)
P <- readRDS(filename_P)
C1 <- readRDS(filename_C1)
C2 <- readRDS(filename_C2)

print(head(T))

print(head(P))

print(head(C1))

print(head(C2))


#print(markers)

if (markers) { C1 <- C1[C2$geneName,]}


#Explicit typecasting for T to run with tibble
T <- data.frame(T)

#Run deconv, params may be exported to config later
res <- CIBERSORT(sig_matrix = C1, mixture_file = T, QN = FALSE, absolute = FALSE, perm = 0)
res <- t(res[,1:(ncol(res)-3)])
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]

#Save and exit
saveRDS(res, file=filename_O)








