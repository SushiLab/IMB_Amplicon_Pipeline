#!/usr/bin/env Rscript
suppressMessages(library(optparse))

# Define arguments
option_list = list(
  make_option(c("-i", "--path_to_seqtab"), type="character", default=NULL,help="Path to the sequence table file (RDS file containing a matrix with sequences as columns and samples as rows)", metavar="character"),
  make_option(c("-s", "--path_to_training_set"), type="character", default=NULL,help="Path to the SILVA training set (will be downloaded if it's not provided)", metavar="character"),
  make_option(c("-c", "--threshold"), type="integer", default=40,help="IdTaxa threshold (default = 40)", metavar="integer"),
  make_option(c("-t", "--threads"), type="integer", default=1,help="Number of threads (default = 1)", metavar="integer"),
  make_option(c("-o", "--out_path"), type="character", default=NULL,help="Path to the output file (table with taxonomy as a tab-delimitted file)", metavar="character")
); 

description<-paste("The program loads an RDS file containing a sequence table and assigns the taxonomy of ASVs/OTUs using IDTAXA\n\n")


opt_parser = OptionParser(option_list=option_list,description = description);
opt = parse_args(opt_parser);

if (is.null(opt$path_to_seqtab) | is.null(opt$out_path)){
  print_help(opt_parser)
  stop("At least one argument must be supplied for -i and -o", call.=FALSE)
}


library(DECIPHER)
library(data.table)
library(tidyverse)

path_to_seqtab<-opt$path_to_seqtab
path_to_training_set<-opt$path_to_training_set
threads<-opt$threads
out_path<-opt$out_path
threshold<-opt$threshold

# Check if the training set exists or download it and load it
if (is.null(path_to_training_set)){
  cat("Training set not provided. It will be downloaded\n")
  system(paste("wget --content-disposition -P ./ http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData",sep=""))
  path_to_training_set<-"SILVA_SSU_r138_2019.RData"
} else{
  cat("Training set already exists. Using local copy\n")
}
load(path_to_training_set)

# Read the RDS file
seqtab<-readRDS(path_to_seqtab)
seqs_fasta<-DNAStringSet(x=as.character(colnames(seqtab)))
names(seqs_fasta)<-as.character(colnames(seqtab))


# Run IDTAXA and parse
annot <- IdTaxa(seqs_fasta, trainingSet=trainingSet, strand="top", processors=threads,threshold=threshold)

annot_df<-sapply(annot,function(x){as.data.frame(x) %>% mutate(annot=paste(rank,taxon,round(confidence,2),sep=";")) %>% summarise(tax=paste(annot,collapse="|"))}) %>%
  unlist() %>%
  as.data.frame() %>%
  rename(tax=".") %>%
  rownames_to_column(var="seq") %>%
  mutate(seq=gsub(".tax$","",seq))

seqtab_annot<-t(seqtab) %>%
  as.data.frame() %>%
  rownames_to_column(var="seq") %>%
  left_join(annot_df,by="seq") %>%
  select(seq,tax,everything())

# Save file
fwrite(seqtab_annot,file=out_path,sep="\t")

