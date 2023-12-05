library(ggplot2)
library(patchwork)
library(dada2)
args = commandArgs(trailingOnly=TRUE)

raw_r1 <- args[1]
raw_r2 <- args[2]
cutadapt_r1 <- args[3]
cutadapt_r2 <- args[4]
filt_r1 <- args[5]
filt_r2 <- args[6]
outfile <- args[7]

# plotQualityProfile plots a visual summary of the distribution of quality scores as a function of sequence position for the input fastq file
topr1 <- plotQualityProfile(raw_r1) + labs(title="R1 - RAW")
topr2 <- plotQualityProfile(raw_r2) + labs(title="R2 - RAW")

midr1 <- plotQualityProfile(cutadapt_r1) + labs(title="R1 - CUTADAPT")
midr2 <- plotQualityProfile(cutadapt_r2) + labs(title="R2 - CUTADAPT")

botr1 <- plotQualityProfile(filt_r1) + labs(title="R1 - QC FILTERED")
botr2 <- plotQualityProfile(filt_r2) + labs(title="R2 - QC FILTERED")



plot <- (topr1 | topr2) /
  (midr1 | midr2) /
  (botr1 | botr2)

ggsave(outfile,plot=plot)
