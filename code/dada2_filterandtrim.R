library(dada2)
packageVersion("dada2")
args = commandArgs(trailingOnly=TRUE)

# paired end
if (length(args) == 12){
    infqgz1 <- args[1] # Input fastq 1
    infqgz2 <- args[2] # Input fastq 2
    outfqgz1 <- args[3] # Output fastq 1
    outfqgz2 <- args[4] # Output fastq 2
    maxee <- as.numeric(args[5]) # After truncation, reads with higher than maxEE "expected errors" will be discarded
    truncq <- as.numeric(args[6]) # Truncate reads at the first instance of a quality score less than or equal to truncQ
    maxn <- as.numeric(args[7]) # Sequences with more than maxN Ns will be discarded
    compress <- as.logical(args[8]) # If output fastq files are gzipped
    minlen <- as.numeric(args[9]) # Remove reads with length less than minLen
    trunclen_r1 <- as.numeric(args[10]) # Truncate reads after truncLen bases. Reads shorter than this are discarded.
    trunclen_r2 <- as.numeric(args[11])
    threads <- as.numeric(args[12])

    # Filter and trim fastq files, outputs trimmed reads which passed the filters
    filterAndTrim(fwd=infqgz1, filt=outfqgz1, rev=infqgz2, filt.rev=outfqgz2, matchIDs=TRUE, maxEE=maxee, truncQ=truncq,
                    maxN=maxn, rm.phix=TRUE, compress=compress, verbose=TRUE, multithread=threads, minLen=minlen,
                    truncLen = c(trunclen_r1, trunclen_r2))
}

# single end
if (length(args) == 9){
    infqgz1 <- args[1] # Input fastq 1
    outfqgz1 <- args[2] # Output fastq 1
    maxee <- as.numeric(args[3]) # After truncation, reads with higher than maxEE "expected errors" will be discarded
    truncq <- as.numeric(args[4]) # Truncate reads at the first instance of a quality score less than or equal to truncQ
    maxn <- as.numeric(args[5]) # Sequences with more than maxN Ns will be discarded
    compress <- as.logical(args[6]) # If output fastq files are gzipped
    minlen <- as.numeric(args[7]) # Remove reads with length less than minLen
    trunclen_r1 <- as.numeric(args[8]) # Truncate reads after truncLen bases. Reads shorter than this are discarded.
    threads <- as.numeric(args[9])

    # Filter and trim fastq files, outputs trimmed reads which passed the filters
    filterAndTrim(fwd=infqgz1, filt=outfqgz1, maxEE=maxee, truncQ=truncq,
                    maxN=maxn, rm.phix=TRUE, compress=compress, verbose=TRUE, multithread=threads, minLen=minlen,
                    truncLen = c(trunclen_r1))
}

