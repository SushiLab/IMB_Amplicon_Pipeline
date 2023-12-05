args = commandArgs(trailingOnly=TRUE)


input.rds <- args[1]

output.stats <- args[2]

seqtab <- readRDS(input.rds)
write.table(rowSums(seqtab), col.names = FALSE, file = output.stats, sep = "\t", quote = FALSE)
