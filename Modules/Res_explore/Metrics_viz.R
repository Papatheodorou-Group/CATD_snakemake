# Metrics visualization script
#
# @zgr2788

#Initialize libs
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(viridis))

sampleName <- commandArgs(trailingOnly = TRUE)

#Grab files
filenames <- list.files("Metrics") %>% lapply(., function(x) { x <- paste0("Metrics/", x) }) %>% lapply(., function(x) { grep(sampleName, x, value = TRUE) }) %>% .[lengths(.)!=0] %>% as.character(.)
files <- lapply(filenames, function(x) readRDS(x) %>% lapply(., function(x) ifelse(is.na(x), 0, x)))
names(files) <- lapply(filenames, function(x) gsub("Metrics/.*res_(.*).rds", "\\1", x))


#Implement visualization for metrics and combinations of metrics below

#Basic Barplot
generateBarPlot <- function(files, metric){
  if (length(files[[metric]]))
  {
    #df <- t(data.frame(files[[metric]]))
    df <- data.frame(as.numeric(files[[metric]]), row.names = names(files[[metric]]))

    # Config y label
    y_lab <- switch(metric,
      "avgcos" = { "Average Cosine Similarity" },
      "rmse" = { "RMSE" },
      "mdist" = { "Minkowski Distance" },
      "dcor" = { "Distance Correlation" },
      "pcor" = { "Pearson Correlation" },
      "mae" = { "MAE" },
      "r2" = { "R2" }
    )


    # Config yticks
    y_axp <- switch(metric,
      "avgcos" = { c(0, 1, 5) },
      "rmse" = { c(0, round(max(df), 1), 5) },
      "mdist" = { c(0, round(max(df), 1), 5) },
      "dcor" = { c(0, 1, 5) },
      "pcor" = { c(0, 1, 5) },
      "mae" = { c(0, round(max(df), 1), 5) },
      "r2" = { c(0, 1, 5) }
    )

    #Get plot
    png(filename = paste0("Plots/",sampleName, "_barplot_", metric, ".png"), width = 1920*3, height = 1080*3, res=300)
    points <- barplot(df[,1], beside=TRUE, col = viridis(nrow(df), direction = -1), yaxp = y_axp)
    mtext(side=1, line=3, "Methods", font=2,cex=1)
    mtext(side=2, line=3, y_lab, font=2,cex=1)
    text(x = points, y = par("usr")[3], labels = rownames(df), xpd = NA, srt = 30, adj = 1, cex = 0.9)
    dev.off()
  }

}

#Generate barplots
for (i in 1:length(files)) { generateBarPlot(files, names(files)[[i]]) }
