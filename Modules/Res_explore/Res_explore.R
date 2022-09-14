#Results exploration script with the given metrics
#
# @zgr2788


#Initialize libs
suppressMessages(library(dplyr))
suppressMessages(library(energy))
suppressMessages(library(lsa))
suppressMessages(library(caret))



#Function to get average cosine similarity over all columns
avg_cossim <- function(P, res) {
  cos <- 0
  for(i in 1:ncol(P)) { cos <- cos + (cosine(P[,i], res[,i])[1]) }
  return (cos/ncol(P))
}


#Get metrics to include in comparison
args <- commandArgs(trailingOnly = TRUE)
filename_P <- args[1]
sampleName <- args[2]
metrics <- args[3:(length(args)-1)]
minkowski_p <- as.numeric(args[length(args)])

P <- readRDS(filename_P)

#Get files as list
filenames <- list.files("Output") %>% lapply(., function(x) { x <- paste0("Output/", x) }) %>% lapply(., function(x) { grep(sampleName, x, value = TRUE) }) %>% .[lengths(.)!=0] %>% as.character(.)
files <- lapply(filenames, readRDS) %>% lapply(., function(x) ifelse(is.na(x), 0, x))
names(files) <- lapply(filenames, function(x) gsub("Output/.*res_(.*).rds", "\\1", x))




#Process metrics
for (i in 1:length(metrics))
{
  metric <- metrics[i]

  switch(metric,

  "rmse"={
    rmse <- lapply(files, function(x) sqrt(mean(as.matrix((x - P)^2))))
    names(rmse) <- names(files)
    rmse <- rmse[order(unlist(rmse),decreasing=FALSE)]
    saveRDS(rmse, paste0("Metrics/", sampleName,"_res_rmse.rds"))
  },

  "mdist"={
    m_dist <- lapply(files, function(x) dist(rbind(as.vector(x), as.vector(unlist(P))), method = "minkowski", p = minkowski_p))
    m_dist <- m_dist[order(unlist(m_dist),decreasing=FALSE)]
    saveRDS(m_dist, paste0("Metrics/", sampleName,"_res_mdist.rds"))
  },

  "dcor"={
    d_cor <- lapply(files, function(x) dcor(P, x))
    d_cor <- d_cor[order(unlist(d_cor),decreasing=TRUE)]
    saveRDS(d_cor, paste0("Metrics/", sampleName,"_res_dcor.rds"))
  },

  "pcor"={
    p_cor <- lapply(files, function(x) cor(c(as.matrix(P)), c(as.matrix(x))))
    p_cor <- p_cor[order(unlist(p_cor),decreasing=TRUE)]
    saveRDS(p_cor, paste0("Metrics/", sampleName,"_res_pcor.rds"))
  },

  "avgcos"={
    cos <- lapply(files, function(x) avg_cossim(P,x))
    cos <- cos[order(unlist(cos),decreasing=TRUE)]
    saveRDS(cos, paste0("Metrics/", sampleName,"_res_avgcos.rds"))
  },

  "mae"={
    mae <- lapply(files, function(x) mean(as.matrix(abs((x - P)))))
    mae <- mae[order(unlist(mae),decreasing=FALSE)]
    saveRDS(mae, paste0("Metrics/", sampleName,"_res_mae.rds"))
  },

  "r2"={
    R2_s <- lapply(files, function(x) R2(c(as.matrix(P)), c(as.matrix(x))))
    R2_s <- R2_s[order(unlist(R2_s),decreasing=TRUE)]
    saveRDS(R2_s, paste0("Metrics/", sampleName,"_res_r2.rds"))
  }
  #Implement more here!
  )
}
