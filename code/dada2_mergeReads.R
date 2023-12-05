library(dada2)
packageVersion("dada2")

# read in arguments
args = commandArgs(trailingOnly=TRUE)
sample_file.r1 <- args[1]
sample_file.r2 <- args[2]
infile.r1 <- args[3]
infile.r2 <- args[4]
outfile.seq_tab.m <- args[5]
outfile.dd.m <- args[6]
orientation <- args[7]

dd.r1 <- readRDS(infile.r1)
dd.r2 <- readRDS(infile.r2)

# validate that orientation is either fw or rev
if (orientation == 'fw') {
  orientation
  } else if (orientation == 'rev') {
    orientation
  } else {
  stop('Orientation can only be fw or rev')
}

# extract sample names
sample.files.r1 <- read.csv(sample_file.r1, header = FALSE, sep = '\t', stringsAsFactors = FALSE)[2]
sample.files.r2 <- read.csv(sample_file.r2, header = FALSE, sep = '\t', stringsAsFactors = FALSE)[2]
s.f.r1 <- sort(sample.files.r1$V2)
s.f.r2 <- sort(sample.files.r2$V2)
sample.names.r1 <- sapply(strsplit(basename(s.f.r1), "_R1"), `[`, 1)
sample.names.r2 <- sapply(strsplit(basename(s.f.r2), "_R2"), `[`, 1)
sample.names.r1
sample.names.r2

# check that file names match
if(!identical(sample.names.r1, sample.names.r2)) stop("Forward and reverse files do not match.")

# assign the names
names(s.f.r1) <- sample.names.r1
names(s.f.r2) <- sample.names.r2

################
## MERGE DATA ##
################
if (orientation == 'fw') {
  mergers <- mergePairs(dd.r1, s.f.r1, dd.r2, s.f.r2, verbose = TRUE)
} else if (orientation == 'rev') {
  mergers <- mergePairs(dd.r2, s.f.r2, dd.r1, s.f.r1, verbose = TRUE)
} else {
  stop('Orientation can only be fw or rev')
}


seq_tab.m <- makeSequenceTable(mergers)
saveRDS(mergers, file = outfile.dd.m)
saveRDS(seq_tab.m, file = outfile.seq_tab.m)
