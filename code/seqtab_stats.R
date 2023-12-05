args = commandArgs(trailingOnly=TRUE)
input.rds <- args[1]
output.stats <- args[2]
seq_tab <- readRDS(input.rds)

# Create stats file from sequence table created by running dada2 inference.
write.table(rowSums(seq_tab), col.names = FALSE, file = output.stats, sep = "\t", quote = FALSE)
