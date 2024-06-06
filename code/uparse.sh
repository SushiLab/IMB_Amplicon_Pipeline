while getopts a:t:o:b:n: flag
do
    case "${flag}" in
        a) input_asvs_f=${OPTARG};;
        t) input_otus_f=${OPTARG};;
        o) output_d=${OPTARG};;
        b) database=${OPTARG};;
        n) project_name=${OPTARG};;
    esac
done

mkdir -p $output_d

# Search for last common ancestor. First, add semicolon to the end of the species name (awk). Secondly, get the last rank that's in common (using sed).
function lca(){ cat $@ | awk '{print $0";"}' |  sed -e '$!{N;s/^\(.*\;\).*\n\1.*$/\1\n\1/;D;}'; return; }

# Perform sequence alignment using usearch
# Search for 500 high identity hits (maxaccepts) in global alignment
usearch -usearch_global $input_asvs_f -db $database -id 0.8 -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out $output_d/"$project_name"_otus.tax -threads 32 &> $output_d/"$project_name"_taxonomy_otus.log
usearch -usearch_global $input_otus_f -db $database -id 0.8 -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out $output_d/"$project_name"_asvs.tax -threads 32 &> $output_d/"$project_name"_taxonomy_asvs.log

# Extract information from the created 'otus.tax' and 'asvs.tax' files and search from last common ancestor ('otus.lca' and 'asvs.lca')
# Get the unique ASVs & OTUs, then get the id (percent identity) and res (the last rank on the tree that's identical between all the hits per ASV/OTU, see example below)

# Example input:
## Bacteria;Firmicutes;Bacilli;Lactobacillales;Listeriaceae;Listeria;Listeria monocytogenes;
## Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecium;
# Example output for the last common ancestor:
## Bacteria;Firmicutes;Bacilli;Lactobacillales;
for i in $(cut -f 1 -d $'\t' $output_d/"$project_name"_otus.tax | sort | uniq); do id=$(grep -m 1 -P $i'\t' $output_d/"$project_name"_otus.tax | cut -f 3 -d$'\t'); res=$(grep -P $i'\t' $output_d/"$project_name"_otus.tax | cut -f 2 -d$'\t' | cut -f 1 -d ' ' --complement | lca); echo -e $i'\t'$res'\t'$id; done > $output_d/"$project_name"_otus.lca
for i in $(cut -f 1 -d $'\t' $output_d/"$project_name"_asvs.tax | sort | uniq); do id=$(grep -m 1 -P $i'\t' $output_d/"$project_name"_asvs.tax | cut -f 3 -d$'\t'); res=$(grep -P $i'\t' $output_d/"$project_name"_asvs.tax | cut -f 2 -d$'\t' | cut -f 1 -d ' ' --complement | lca); echo -e $i'\t'$res'\t'$id; done > $output_d/"$project_name"_asvs.lca

