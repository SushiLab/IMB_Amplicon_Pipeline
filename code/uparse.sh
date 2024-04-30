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

# Search for last common ancestor
function lca(){ cat $@ | sed -e '$!{N;s/^\(.*\).*\n\1.*$/\1\n\1/;D;}' | awk -F ";" '{$NF=""; OFS=";"; print $0}'; return; }

# Perform sequence alignment using usearch
# Search for 500 high identity hits (maxaccepts) in global alignment
usearch -usearch_global $input_asvs_f -db $database -id 0.8 -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out $output_d/"$project_name"_otus.tax -threads 32 &> $output_d/"$project_name"_taxonomy_otus.log
usearch -usearch_global $input_otus_f -db $database -id 0.8 -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out $output_d/"$project_name"_asvs.tax -threads 32 &> $output_d/"$project_name"_taxonomy_asvs.log

# Extract information from the created 'otus.tax' and 'asvs.tax' files and search from last common ancestor ('otus.lca' and 'asvs.lca')
for i in $(cut -f 1 -d $'\t' $output_d/"$project_name"_otus.tax | sort | uniq); do id=$(grep -m 1 -P $i'\t' $output_d/"$project_name"_otus.tax | cut -f 3 -d$'\t'); res=$(grep -P $i'\t' $output_d/"$project_name"_otus.tax | cut -f 2 -d$'\t' | cut -f 1 -d ' ' --complement | lca); echo -e $i'\t'$res'\t'$id; done > $output_d/"$project_name"_otus.lca
for i in $(cut -f 1 -d $'\t' $output_d/"$project_name"_asvs.tax | sort | uniq); do id=$(grep -m 1 -P $i'\t' $output_d/"$project_name"_asvs.tax | cut -f 3 -d$'\t'); res=$(grep -P $i'\t' $output_d/"$project_name"_asvs.tax | cut -f 2 -d$'\t' | cut -f 1 -d ' ' --complement | lca); echo -e $i'\t'$res'\t'$id; done > $output_d/"$project_name"_asvs.lca

