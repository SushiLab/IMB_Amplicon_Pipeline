library(dada2)
packageVersion("dada2")
library(readr)

# Read in arguments
args = commandArgs(trailingOnly=TRUE)
sample_file <- args[1]
err.rds <- args[2]
outfile.tab <- args[3]
outfile.dd <- args[4]
threads <- as.numeric(args[5])
ref_sequence_path <- args[6]
script_folder <- args[7]

# Read in dada2_learnErrors.R to use the function loess_error_function_mod
source(paste0(script_folder, "loess_error_function.R"))

# Read in sample file
sample.table <- read.csv(sample_file, header=FALSE, sep='\t', stringsAsFactors = FALSE)
sample.names <- sample.table$V1
sample.files <- sample.table$V2
s.f <- sample.files
names(s.f) <- sample.names
s.f <- sort(s.f)

# Read in error rate from file created by dada2_learnErrors
err <- readRDS(err.rds)

# Perform inference
if(ref_sequence_path != 'None'){
    # Depending on the presence of a reference sequence file, read in sequences and perform DADA2 analysis.
    lines = read_lines(ref_sequence_path)
    header_index = which(grepl(">", lines))
    sequences = unlist(lapply(split(lines, rep(header_index, diff(c(header_index, length(lines) + 1)))), function(x) paste(x[2:length(x)], collapse="")))
    names(sequences) <- sapply(lines[header_index], function(x) sub(">", "", x))
    # Run DADA2 with prior information and pseudo-pooling
    dd <- dada(s.f, err = err, pool = 'pseudo', multithread = threads, errorEstimationFunction = loess_error_function_mod, priors = sequences)
}else{
    dd <- dada(s.f, err = err, pool = 'pseudo', multithread = threads, errorEstimationFunction = loess_error_function_mod)
}

# Save
dd[[1]]
# Generate sequence table from DADA2 object
seqtab <- makeSequenceTable(dd)
saveRDS(seqtab, file = outfile.tab)
saveRDS(dd, file = outfile.dd)

