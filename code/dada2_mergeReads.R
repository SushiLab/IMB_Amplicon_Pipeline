library(dada2); packageVersion("dada2")
args = commandArgs(trailingOnly=TRUE)

#args <- c('/Users/hans/Desktop/fw_R1.learnerrors.samples', 
#          '/Users/hans/Desktop/fw_R2.learnerrors.samples',
#          '/Volumes/biol_micro_sunagawa/Projects/PAN/TARA_PACIFIC_METAB_PAN/code/pipeline/fwout.r1.rds','/Volumes/biol_micro_sunagawa/Projects/PAN/TARA_PACIFIC_METAB_PAN/code/pipeline/fwout.r2.rds', '/Volumes/biol_micro_sunagawa/Projects/PAN/TARA_PACIFIC_METAB_PAN/code/pipeline/fwout.m.rds', 
#          '4', 'fw')

samplefile.r1 <- args[1]
samplefile.r2 <- args[2]


infile.r1 <- args[3]
infile.r2 <- args[4]

outfile.seqtab.m <- args[5]
outfile.dd.m <- args[6]

orientation <- args[7]

if (orientation == 'fw') {
  orientation
  } else if (orientation == 'rev') {
    orientation
  } else {
  stop('Orientation can only be fw or rev')
}


sample.files.r1 <- read.csv(samplefile.r1, header=FALSE, sep='\t', stringsAsFactors = FALSE)[2]
sample.files.r2 <- read.csv(samplefile.r2, header=FALSE, sep='\t', stringsAsFactors = FALSE)[2]

s.f.r1 <- sort(sample.files.r1$V2)
s.f.r2 <- sort(sample.files.r2$V2)


sample.names.r1 <- sapply(strsplit(basename(s.f.r1), "_R1"), `[`, 1)
sample.names.r2 <- sapply(strsplit(basename(s.f.r2), "_R2"), `[`, 1)

sample.names.r1
sample.names.r2
if(!identical(sample.names.r1, sample.names.r2)) stop("Forward and reverse files do not match.")


names(s.f.r1) <- sample.names.r1
names(s.f.r2) <- sample.names.r2


dd.r1 <- readRDS(infile.r1)
dd.r2 <- readRDS(infile.r2)




if (orientation == 'fw') {
  mergers <- mergePairs(dd.r1, s.f.r1, dd.r2, s.f.r2, verbose = TRUE)
} else if (orientation == 'rev') {
  mergers <- mergePairs(dd.r2, s.f.r2, dd.r1, s.f.r1, verbose = TRUE)
} else {
  stop('Orientation can only be fw or rev')
}


seqtab.m <- makeSequenceTable(mergers)
saveRDS(mergers, file = outfile.dd.m)
saveRDS(seqtab.m, file = outfile.seqtab.m)
