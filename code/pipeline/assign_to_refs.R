cat("Testing USEARCH is accessible...")
usearch_test <- try(expr = system(command = "usearch",intern = T),silent = T)
if(class(usearch_test) == "try-error"){
      stop("USEARCH executable is not found. It should be loaded and be executable with the exact command `'usearch'. You can add aliases to your '.bashrc' file.", call.=FALSE)
}
cat("DONE\n")

args = commandArgs(trailingOnly=TRUE)

#args <- c("scratch/BM10_library/BM10_library.asvs.tsv", "scratch/BM10_library/BM10_library.asvs.fasta", "lib/atlsphere_16s_trimmed.fasta", 4, "scratch/BM10_library/BM10_library.refs.tsv")

asvtab_path <- args[1]
asvs_fasta <- args[2]
refs_fasta <- args[3]
alntab_path <- sub(".fasta", ".aln", asvs_fasta)
maxe <- as.numeric(args[4])
output_path <- args[5]

cat("Running USEARCH algorithm...")
cmd <- system(command = paste("usearch -usearch_global", asvs_fasta, "-db",  refs_fasta, "-strand plus -blast6out", alntab_path, "-id 0.0 -maxaccepts 0 -threads 32"))
cat("DONE\n")

asvtab <- read.table(asvtab_path, sep="\t", header=T)
alntab <- read.table(alntab_path, sep="\t")
alntab <- alntab[alntab[,5] <= maxe,]

assign <- function(alntab){
    alntab <- alntab[order(alntab[,3], decreasing=T),]
    best_hits <- which(alntab[,3] == max(alntab[,3]))
    if(length(best_hits) > 1){
        best_hits <- which(alntab[,4] == max(alntab[,4]))
        if(length(best_hits) > 1){
            return("ambiguous")
        }else{
            return(alntab[best_hits,2])
        }
    }else{
        return(alntab[best_hits,2])
    }
}

assignments <- unlist(lapply(split(alntab, alntab[,1]), assign))

asvs <- unlist(lapply(asvtab[,1], function(x) strsplit(x, ";")[[1]][1]))
asvtab$assignment <- assignments[asvtab[,1]]
asvtab <- asvtab[,c(1:5, ncol(asvtab), 6:(ncol(asvtab)-1))]
found_refs <- unique(asvtab$assignment[!is.na(asvtab$assignment)])
found_refs <- found_refs[found_refs != "ambiguous"]
asvtab$assignment[is.na(asvtab$assignment)] <- "none"

# collapse identical assignments
crows <- lapply(found_refs, function(x) colSums(asvtab[asvtab$assignment == x, 7:ncol(asvtab)]))
reftab <- do.call(rbind, crows)
rownames(reftab) <- found_refs
ambitab <- asvtab[asvtab$assignment == "ambiguous", 7:ncol(asvtab)]
rownames(ambitab) <- paste(asvs[asvtab$assignment == "ambiguous"], "(ambiguous)")
nonetab <- asvtab[asvtab$assignment == "none", 7:ncol(asvtab)]
rownames(nonetab) <- paste(asvs[asvtab$assignment == "none"], "(none)")
fintab <- rbind(reftab, ambitab, nonetab)

write.table(fintab, output_path, sep="\t", quote=F)
