library(dada2); packageVersion("dada2")

# Define and parse parameters
args = commandArgs(trailingOnly=TRUE)

wbim.file <- args[1]
nobim.file <- args[2]
threads <- as.numeric(args[3])

wbim.tab <- readRDS(wbim.file)

nobim.tab <- removeBimeraDenovo(wbim.tab, method="pooled", multithread=threads, verbose=TRUE)

saveRDS(nobim.tab, file = nobim.file)



