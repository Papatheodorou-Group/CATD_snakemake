# script to get and plot average of DWLS, FARDEEP and EpiDISH
# and calculate average correlation between the three methods
# @nadjano

suppressMessages(library(dplyr))
#install.packages("Polychrome", repos='http://cran.us.r-project.org')
library(Polychrome)

CiberBarFrazer <- function(ciber, colors, main0, legy){
    ## function to get barplot per sample
    ## from https://github.com/mkrdonovan/gtex_deconvolution
    par(mar = c(2, 4, 2, 0.5))
    
    nsamples =nrow(ciber)
    ciber = as.data.frame(t(ciber) * 100)
    
    ciber$color = colors[ match(colors$celltypes, rownames(ciber)),  "color"]
    barplot(as.matrix(ciber[, seq(1, (ncol(ciber) - 1))]), 
            las = 2, 
            col = ciber$color, 
            border=NA, 
            names.arg = rep("", ncol(ciber) - 1), 
            ylab = "Fraction clusters", 
            main = main0, 
            space=0, 
            cex.main = 1)

    text(nsamples * .05, 90, paste("n = ", nsamples, sep = ""), pos = 4, cex = 1)

    legend(nsamples * .05, legy, gsub("_", " ", colors$celltypes), bty = "n",
           pch = rep(22, nrow(colors)),
           pt.cex = rep(4, nrow(colors)),
           pt.bg = colors$color,
           y.intersp = 1, cex = 1.1
          )   
}
sampleName <- commandArgs(trailingOnly = TRUE)

#Get files as list 
#only select DWLS, EpiDISH and FARDEEP proportions for consensus
filenames <- list.files("Output/", pattern="DWLS|EpiDISH|FARDEEP") 
            %>% lapply(., function(x) { x <- paste0("Output/", x) }) 
            %>% lapply(., function(x) { grep(sampleName, x, value = TRUE) }) 
            %>% .[lengths(.)!=0] %>% as.character(.)

files <- lapply(filenames, readRDS)
names(files) <- lapply(filenames, function(x) gsub("Output/.*res_(.*).rds", "\\1", x))

for (i in 1:length(files)){
    files[[i]] <- files[[i]][,order(colnames(files[[i]]))]
    files[[i]] <- files[[i]][order(rownames(files[[i]])),]
  }
  # Calculate Pearson correlation for each RUN between the three methods
  # to see if the difference between methods is too big
  mean_vector <- numeric()
  for (i in 1:dim(files[[1]])[2]) {
    cor_mat <- cor(data.frame(files[[1]][, i], files[[2]][, i], files[[3]][, i]))
    
    # Set diagonal to NA
    diag(cor_mat) <- NA
    
    mean_value <- mean(cor_mat, na.rm = TRUE)
    mean_vector[i] <- mean_value
  }
  
  mean_corr_all = round(mean(mean_vector), 3)
  # Check mean correlation
  if (mean(mean_vector) > 0.6) {
    cat(paste0('mean_correlation:', round(mean(mean_vector), 3)))
  } else {
    cat(paste0('correlation_between_methods_lower_than_0.6'))
  }

#get average proportions per sample and celltype accros the three methods
sum_props = Reduce('+',files)
prop = sum_props/length(files)
#Save consensus proportions
saveRDS(prop, paste0("Consensus/",sampleName, "_consensus.rds"))


#get sd
vec <- unlist(files, use.names = TRUE)
DIM <- dim(files[[1]])
n <- length(files)

list.mean <- tapply(vec, rep(1:prod(DIM),times = n), mean)
attr(list.mean, "dim") <- DIM
list.mean <- as.data.frame(list.mean)

list.sd <- tapply(vec, rep(1:prod(DIM),times = n), sd)
attr(list.sd, "dim") <- DIM
list.sd <- as.data.frame(list.sd)
colnames(list.sd) = colnames(files[[1]])

sd_norm = rowMeans(list.sd)/rowMeans(list.mean)
prop = data.frame(t(prop))

#add sd to columnames 
colnames(prop) = paste0(colnames(prop), '(sd=', 100* round(rowMeans(list.sd),3), ')')
#order to have largest proportions first
top = names(prop[order(colSums(prop), decreasing = T)][1])
second_top = names(prop[order(colSums(prop), decreasing = T)][2])
third_top = names(prop[order(colSums(prop), decreasing = T)][3])
color = data.frame(celltypes = colnames(prop),
                   name      =  colnames(prop),
                   color     =  sky.colors(ncol(prop)))

png(filename = paste0("ConsensusPlot/",sampleName, "_consensus.png"), 
    width = 500*3, 
    height = 500*3, 
    res=300)
CiberBarFrazer(prop[order(-prop[top], -prop[second_top], -prop[third_top]),], color, paste(sampleName, '_', mean_corr_all), 90)
dev.off()
