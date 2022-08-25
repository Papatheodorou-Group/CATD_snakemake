# Script for violinplots for cosine sim and rmse
#
# @zgr2788

#libs
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(lsa))




#Get metrics
args <- commandArgs(trailingOnly = TRUE)
filename_P <- args[1]
sampleName <- args[2]
metrics <- args[3:(length(args)-1)]
minkowski_p <- as.numeric(args[length(args)])

P <- readRDS(filename_P)

#Get files as list
filenames <- list.files("Output") %>% lapply(., function(x) { x <- paste0("Output/", x) }) %>% lapply(., function(x) { grep(sampleName, x, value = TRUE) }) %>% .[lengths(.)!=0] %>% as.character(.)
files <- lapply(filenames, readRDS)
names(files) <- lapply(filenames, function(x) gsub("Output/.*res_(.*).rds", "\\1", x))

getPlot <- function(df, name) {
  ggplot(df, aes(x = Methods, y = Value, fill = Methods)) + geom_violin(width = 1.3) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0("Plots/", sampleName, "_vioplot_", name, ".png"), width = 1920*3, height = 1080*3, dpi=300, units = "px")
}


#Process metrics
for (i in 1:length(metrics))
{
  metric <- metrics[i]

  switch(metric,

  "rmse"={
    rmse <- lapply(files, function(x) sqrt(colSums(as.matrix((x - P)^2)) / nrow(P)))
    rmse <- data.frame(rmse)
    rmse <- rmse %>% gather(Methods, Value)
    getPlot(rmse, metric)
  },

  "avgcos"={
    cos <- lapply(files, function (x)  { mapply(function(y, z) cosine(y,z), data.frame(x), P) })
    cos <- data.frame(cos)
    cos <- cos %>% gather(Methods, Value)
    getPlot(cos, metric)
  }
  )
}
