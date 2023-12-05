library(dada2)
packageVersion("dada2")
library(ggplot2)

# read in arguments
args = commandArgs(trailingOnly=TRUE)
sample_file <- args[1]
out_file <- args[2]
threads <- as.numeric(args[3])
nbases <- as.numeric(args[4])
script_folder <- args[5]

# read in the loess error function
source(paste0(script_folder, "loess_error_function.R"))

# read in sample files
sample.files <- read.csv(sample_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)[2]
s.f <- sample.files$V2

# create path for out_file
out_file.plot <- paste(out_file, '.pdf', sep = '')
out_file.plot

# learn the errors
err <- learnErrors(s.f, nbases = nbases, multithread = threads, randomize = TRUE, verbose = 1,
                    errorEstimationFunction = loess_error_function_mod)
saveRDS(err, file = out_file)

# plot the errors
plot <- plotErrors(err, nominalQ = TRUE)
ggsave(out_file.plot, plot = plot)
