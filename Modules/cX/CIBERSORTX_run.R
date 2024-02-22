# convert_cibersort_output.R



args <- commandArgs(trailingOnly = TRUE)

filename_T <- args[1]
filename_P <- args[2]
filename_C0 <- args[3]
filename_phenData <- args[4]
filename_input <- args[5]
filename_0 <- args[6]

P <- readRDS(filename_P)

# Read CIBERSORTX output from text file
cibersort_output <- read.table(filename_input, header = TRUE, sep = "\t", row.names = 1)
cibersort_output<- cibersort_output[1:(length(cibersort_output)-3)]
# Convert to a data frame or any necessary processing
res<- t(cibersort_output)
#Save the data frame as an RDS file 
if (filename_P != 'Modules/Psuedobulk/dummy_props.rds') res <- res[order(match(rownames(res), rownames(P))),]



saveRDS(res, file = filename_0)
