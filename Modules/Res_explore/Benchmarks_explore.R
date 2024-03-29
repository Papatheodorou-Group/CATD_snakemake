# Exploration of different modules' benchmarks
#
# @zgr2788


#Libs go here
suppressMessages(library(dplyr))
suppressMessages(library(viridis))


sampleName <- commandArgs(trailingOnly = TRUE)

#Get benchmarks as list
filenames <- list.files("Benchmarks") %>% lapply(., function(x) { x <- paste0("Benchmarks/", x) }) %>% lapply(., function(x) { grep(sampleName, x, value = TRUE) }) %>% .[lengths(.)!=0] %>% as.character(.)
files <- lapply(filenames, function(x) read.csv(x, sep = "\t"))
names(files) <- lapply(filenames, function(x) gsub("Benchmarks/.*_(.*)_benchmark.txt", "\\1", x))
files <- lapply(files, function(x) x <<- x[,c(1,3,9,10)])

#Build dataframe for exploration
benchDf <- data.frame(row.names=names(files[[1]]))
for (i in 1:length(files)) { benchDf[,names(files[i])] <- t(files[[i]]) }
                
                
# Set up the filename for SVG output
output_filename <- paste0(sampleName, "_benchmarks_summarized.pdf")

# Use svg() for vectorized output
pdf(file = output_filename)                
                
                
                
#Get plot
points <- barplot(t(sqrt(as.matrix(benchDf))), beside = TRUE, ylab = "value (sqrt transformed)", main = "Benchmark Parameters", col = viridis(length(files)), yaxp=c(round(min(sqrt(benchDf))),round(max(sqrt(benchDf))),20), horiz = FALSE,width=0.8)
shiftLab <- unname(sapply(c(replicate(nrow(benchDf), colnames(benchDf))), function(x) nchar(x)))
text(x = points, y = c(t(sqrt(as.matrix(benchDf)))) + (1.65*(shiftLab + 5/shiftLab)), labels = c(replicate(nrow(benchDf), colnames(benchDf))), xpd = NA, srt = 90, adj = 1, cex = 5)

#Save plot
dev.off()

                          
