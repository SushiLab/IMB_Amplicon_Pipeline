library(dada2); packageVersion("dada2")
library(readr)
args = commandArgs(trailingOnly=TRUE)

#args <- c('/Users/hans/Desktop/fw_R1.learnerrors.samples', 
#          '/Volumes/biol_micro_sunagawa/Projects/PAN/TARA_PACIFIC_METAB_PAN/data/processed/pipeline_test/pipeline/3learnErrors/H32KCBCX2/1/fw_R1.learnerrors.rds', 
#          '/Volumes/biol_micro_sunagawa/Projects/PAN/TARA_PACIFIC_METAB_PAN/code/pipeline/seqtab.r1.rds', 
#          ''/Volumes/biol_micro_sunagawa/Projects/PAN/TARA_PACIFIC_METAB_PAN/code/pipeline/dd.r1.rds''
#          '4')

samplefile <- args[1]

err.rds <- args[2]

outfile.tab <- args[3]
outfile.dd <- args[4]
threads <- as.numeric(args[5])
ref_sequence_path <- args[6]



loessErrfun_mod <- function (trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A",
                                                "C", "G", "T")), ])
        rlogp <- log10((errs + 1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot,
                         rlogp = rlogp)
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-07
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ], est[4,], 1 - colSums(est[4:6, ]), est[5:6, ], est[7:8, ], 1 -
                   colSums(est[7:9, ]), est[9, ], est[10:12, ], 1 - colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4),
                          "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  return(err)
}













#sample.files <- read.csv(samplefile, header=FALSE, sep='\t', stringsAsFactors = FALSE)[2]
sample.table <- read.csv(samplefile, header=FALSE, sep='\t', stringsAsFactors = FALSE)
sample.files <- sample.table$V2
#s.f <- sort(sample.files)
s.f <- sample.files

#sample.names <- sapply(strsplit(basename(s.f), "_R"), `[`, 1)
sample.names <- sample.table$V1

#if(!identical(sample.names.r1, sample.names.r2)) stop("Forward and reverse files do not match.")


names(s.f) <- sample.names
s.f <- sort(s.f)


err <- readRDS(err.rds)


if(ref_sequence_path != 'None'){
    # Read in sequences
    lines = read_lines(ref_sequence_path)
    header_index = which(grepl(">", lines))
    sequences = unlist(lapply(split(lines, rep(header_index, diff(c(header_index, length(lines)+1)))), function(x) paste(x[2:length(x)], collapse="")))
    names(sequences) <- sapply(lines[header_index], function(x) sub(">", "", x))
    dd <- dada(s.f, err=err, pool='pseudo', multithread = threads, errorEstimationFunction = loessErrfun_mod, priors = sequences)
}else{
    dd <- dada(s.f, err=err, pool='pseudo', multithread = threads, errorEstimationFunction = loessErrfun_mod)
}

dd[[1]]
seqtab <- makeSequenceTable(dd)
saveRDS(seqtab, file = outfile.tab)
saveRDS(dd, file = outfile.dd)

