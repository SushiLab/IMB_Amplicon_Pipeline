library(dada2)
packageVersion("dada2")
library(ggplot2)

# Read in arguments
args = commandArgs(trailingOnly=TRUE)
sample_file <- args[1]
out_file <- args[2]
threads <- as.numeric(args[3])
n_bases <- as.numeric(args[4])
script_folder <- args[5]

# Read in the loess error function
source(paste0(script_folder, "loess_error_function.R"))

# Read in sample files
sample_files <- read.csv(sample_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)[2]
s_f <- sample_files$V2

# Create path for out_file
out_file.plot <- paste(out_file, '.pdf', sep = '')
out_file.plot

# Learn the errors
err <- learnErrors(s_f, nbases = n_bases, multithread = threads, randomize = TRUE, verbose = 1,
                    errorEstimationFunction = loess_error_function_mod)
saveRDS(err, file = out_file)

# Plot the errors
plot <- plotErrors(err, nominalQ = TRUE)
ggsave(out_file.plot, plot = plot)
