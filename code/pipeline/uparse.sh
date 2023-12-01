input_d="scratch"
output_d="scratch/8uparsetax"
db =

mkdir -p $output_d

function lca(){ cat $@ | sed -e '$!{N;s/^\(.*\).*\n\1.*$/\1\n\1/;D;}' | awk -F ";" '{$NF=""; OFS=";"; print $0}'; return; }

usearch -usearch_global $input_d/scratch.otus.fasta -db /nfs/cds/Databases/SILVA/SILVA138/SILVA_138.1_SSURef_NR99_tax_silva.fasta -id 0.8 -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out $output_d/otus.tax -threads 32 &> $output_d/taxonomy.log

for i in $(cut -f 1 -d $'\t' $output_d/otus.tax | sort | uniq); do id=$(grep -m 1 -P $i'\t' $output_d/otus.tax | cut -f 3 -d$'\t'); res=$(grep -P $i'\t' $output_d/otus.tax | cut -f 2 -d$'\t' | cut -f 1 -d ' ' --complement | lca); echo -e $i'\t'$res'\t'$id; done > $output_d/otus.lca

#usearch -otutab $output_d/pool.primermatch.fasta -otus $output_d/otus_uparse.fa -strand both -id 0.97 -otutabout $output_d/otutab.txt -threads 16 &> $output_d/otutab.log

#function lca(){ cat $@ | sed -e '$!{N;s/^\(.*\).*\n\1.*$/\1\n\1/;D;}' | awk -F ";" '{$NF=""; OFS=";"; print $0}'; return; }
#usearch -usearch_global {input.input_d}/scratch.otus.fasta -db /nfs/cds/Databases/SILVA/SILVA138/SILVA_138.1_SSURef_NR99_tax_silva.fasta -id 0.8 -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out {output.output_d}/otus.tax -threads 32 &> {output.output_d}/taxonomy.log
#for i in $(cut -f 1 -d $"\t" {output.output_d}/otus.tax | sort | uniq); do id=$(grep -m 1 -P $i"\t" {output.output_d}/otus.tax | cut -f 3 -d$"\t"); res=$(grep -P $i"\t" {output.output_d}/otus.tax | cut -f 2 -d$"\t" | cut -f 1 -d ' ' --complement | lca); echo -e $i"\t"$res"\t"$id; done > {output.output_d}/otus.lca