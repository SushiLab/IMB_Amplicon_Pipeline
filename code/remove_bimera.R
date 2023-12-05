library(dada2)
packageVersion("dada2")

# read in arguments
args = commandArgs(trailingOnly=TRUE)
wbim.file <- args[1]
nobim.file <- args[2]
threads <- as.numeric(args[3])

wbim.tab <- readRDS(wbim.file)

# Remove bimeric sequences that are formed de novo during the sequencing process
nobim.tab <- removeBimeraDenovo(wbim.tab, method="pooled", multithread=threads, verbose=TRUE)

saveRDS(nobim.tab, file = nobim.file)



