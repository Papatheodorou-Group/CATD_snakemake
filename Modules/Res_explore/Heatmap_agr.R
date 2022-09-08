# Heatmap for method agreements
#
# @zgr2788

suppressMessages(library(dplyr))
suppressMessages(library(lsa))
suppressMessages(library(viridis))



sampleName <- commandArgs(trailingOnly = TRUE)


function(x,y){
  idx <- seq_len(ncol(x))
  distances <- sapply(idx, function(z) cosine(x[,z], y[,z]))
  return(mean(distances))
}





#Get files as list
filenames <- list.files("Output") %>% lapply(., function(x) { x <- paste0("Output/", x) }) %>% lapply(., function(x) { grep(sampleName, x, value = TRUE) }) %>% .[lengths(.)!=0] %>% as.character(.)
files <- lapply(filenames, readRDS)
names(files) <- lapply(filenames, function(x) gsub("Output/.*res_(.*).rds", "\\1", x))

#Reorder files in case this was without proportions
files <- lapply(files, function(x) x[order(match(rownames(x), rownames(files[[1]]))),])

#Generate combs to iterate over
combs <- combn(files, 2)
colnames(combs) <- combn(names(files), 2, paste0, collapse="@")
avgDists <- sapply(seq_len(ncol(combs)), function(z)  as.numeric(table(diag(cor(combs[,z][[1]], combs[,z][[2]])) > 0.9)["TRUE"]) )
#avgDists <- sapply(seq_len(ncol(combs)), function(z)  getAvgDist(combs[,z][[1]], combs[,z][[2]]) )
#avgDists <- sapply(seq_len(ncol(combs)), function(z) mean(diag(cor(combs[,z][[1]], combs[,z][[2]]))) )
names(avgDists) <- colnames(combs)
avgDists[is.na(avgDists)] <- 0


#Generate matrix
mat <- matrix(0, nrow = length(files), ncol = length(files))
diag(mat) <- ncol(combs[,1][[1]])
colnames(mat) <- names(files)
rownames(mat) <- names(files)

#Copy upper triangle to lower
mat[lower.tri(mat, diag=FALSE)] <- avgDists
mat <- t(mat)
mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

#Save plot
png(filename = paste0("Plots/",sampleName, "_metricsHeatmap.png"), width = 1920*3, height = 1080*3, res=300)
heatmap(mat, col = rev(viridis(256)), scale = "none")
dev.off()
