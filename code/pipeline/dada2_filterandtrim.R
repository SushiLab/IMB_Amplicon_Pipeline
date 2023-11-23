library(dada2);
packageVersion("dada2")
args = commandArgs(trailingOnly=TRUE)

infqgz1 <- args[1]
infqgz2 <- args[2]
outfqgz1 <- args[3]
outfqgz2 <- args[4]
maxee <- as.numeric(args[5])
truncq <- as.numeric(args[6])
maxn <- as.numeric(args[7])
compress <- as.logical(args[8])
minlen <- as.numeric(args[9])
trunclen_r1 <- as.numeric(args[10])
trunclen_r2 <- as.numeric(args[11])
threads <- as.numeric(args[12])


filterAndTrim(fwd=infqgz1, filt=outfqgz1, rev=infqgz2, filt.rev=outfqgz2, matchIDs=TRUE, maxEE=maxee, truncQ=truncq, maxN=maxn, rm.phix=TRUE, compress=compress, verbose=TRUE, multithread=threads, minLen=minlen, truncLen = c(trunclen_r1, trunclen_r2))
