#!/usr/bin/env Rscript
suppressMessages(library(optparse))

# Define arguments
option_list = list(
  make_option(c("-i", "--path_to_table"), type="character", default=NULL,help="Path to the ASV table", metavar="character"),
  make_option(c("-f", "--out_path_fasta"), type="character", default="asvs.fasta",help="Path to the FASTA output file", metavar="character"),
  make_option(c("-a", "--out_path_asv"), type="character", default="asvtab.tsv",help="Path to the ASV output file", metavar="character"),
  make_option(c("-o", "--out_path_otu"), type="character", default="otutab.tsv",help="Path to the OTU output file", metavar="character"),
  make_option(c("-u", "--out_path_uparse"), type="character", default="uparse.tsv",help="Path to the uprase file", metavar="character"),
  make_option(c("-x", "--out_path_otu_fasta"), type="character", default="otus.fasta",help="Path to the OTU fasta file", metavar="character")
); 

description<-paste("The program loads an ASV table and uses the UPARSE algorithm to cluster ASVs into OTUs and produce an OTU table.\n\n",
                   "The ASV table (specified in -i) needs to be a tab-delimited file with ASVs as rows and samples as columns. It also needs to have a mandatory column named as 'seq' containing the ASV sequence. All numeric variables should correspond to and only to the abundances in each sample. Any other non-numeric variable is allowed and will be recycled into the OTU table.\n",
                   "From the ASV sequence a FASTA file (specified in -f) is produced and subsequently used to run the UPARSE algorithm. A new ASV table including the UPARSE information is saved (specified in -a). An OTU table is produced and saved (specified in -o). For each OTU, information  from the representative ASV is kept. ASVs labelled as chimeric sequences by UPARSE are kept in the new ASV table but excluded in the OTU table.\n",
                   "USEARCH needs to be accessible and be executable with the exact command 'usearch'.\n")


opt_parser = OptionParser(option_list=option_list,description = description);
opt = parse_args(opt_parser);

if (is.null(opt$path_to_table)){
  print_help(opt_parser)
  stop("At least one argument must be supplied for -i", call.=FALSE)
}

cat("Testing USEARCH is accessible...")
usearch_test<-try(expr = system(command = "usearch",intern = T),silent = T)
if(class(usearch_test)=="try-error"){
  stop("USEARCH executable is not found. It should be loaded and be executable with the exact command `'usearch'. You can add aliases to your '.bashrc' file.", call.=FALSE)
}
cat("DONE\n")

cat("Loading libraries...")
library(data.table)
library(tidyverse)
library(Biostrings)
cat("DONE\n")

path_to_table<-opt$path_to_table
out_path_fasta<-opt$out_path_fasta
out_path_asv<-opt$out_path_asv
out_path_otu<-opt$out_path_otu
out_path_otu_fasta<-opt$out_path_otu_fasta
out_path_uparse<-opt$out_path_uparse

# Load the ASV table
cat("Reading ASV table...")
if (file.exists(path_to_table)){
  asvtab<-fread(path_to_table,sep="\t",header=T,data.table = F)
} else {
  cat("Error: ASV table file does not exist\n")
}
cat("DONE\n")


# Create ASV ID
cat("Creating ASV identifiers...")
asvtab<-asvtab %>% # create ASV id with size annotation
  mutate(total_ab = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(total_ab)) %>%
  mutate(asv=paste("asv_",str_pad(seq(1,nrow(asvtab)),nchar(nrow(asvtab)),pad="0"),";size=",total_ab,sep="")) %>%
  select(-total_ab) %>%
  select(asv,everything())
cat("DONE\n")

# Produce the ASV FASTA file
cat("Generating ASV FASTA file...")
seqs<-DNAStringSet(asvtab$seq)
names(seqs)<-asvtab$asv

# Save the ASV FASTA file
writeXStringSet(seqs,filepath = out_path_fasta)
cat("DONE\n")

# Run UPARSE algorithm
cat("Running UPARSE algorithm...")
cmmd<-system(command = paste("usearch -cluster_otus ",out_path_fasta," -minsize 1 -otus ",out_path_otu_fasta," -uparseout ",out_path_uparse," -relabel otu",sep=""))
cat("DONE\n")

# Parse UPARSE output and add it to the ASV table
cat("Parsing UPARSE output...")
uparseout<-fread(out_path_uparse,sep="\t",header=F,data.table = F) %>%
  separate(V3,into=c("dqt","top","dqm","div","segs","parents",NA),sep=";",remove=F,fill="right") %>%
  mutate(top=gsub("top=","",top)) %>%
  separate(top,into=c("top",NA),sep="\\(",fill="right") %>%
  mutate(V2=ifelse(V2 =="match",top,V2)) %>%
  mutate(dqt=gsub("top=","",dqt)) %>%
  separate(dqt,into=c("dqt",NA),sep="\\(",fill="right") %>%
  mutate(V2=ifelse(V2 =="perfect",dqt,V2)) %>%
  select(asv=V1,otu=V2,uparse_info=V3)
cat("DONE\n")

# Save the ASV table with ASV names and UPARSE info
cat("Saving ASV table with UPARSE info...")
asvtab<-asvtab %>%
  left_join(uparseout,by="asv") %>%
  select(asv,otu,uparse_info,everything())

fwrite(asvtab,file = out_path_asv,sep="\t")
cat("DONE\n")

# Produce and save the OTU table
cat("Generating OTU table...")
otutab<-asvtab %>%
  filter(grepl("^otu",otu)) %>%
  group_by(otu) %>%
  summarise(across(where(is.numeric),sum))

tmp<-asvtab %>%
  select_if(negate(is.numeric)) %>%
  filter(grepl("^otu",otu)) %>%
  select(-asv,-uparse_info) %>%
  group_by(otu) %>%
  slice_head(n = 1)

otutab<-otutab %>%
  left_join(tmp,by="otu") %>%
  select(colnames(tmp),everything())
            
fwrite(otutab,file = out_path_otu,sep="\t")
cat("DONE\n")
cat("FINISHED\n")


