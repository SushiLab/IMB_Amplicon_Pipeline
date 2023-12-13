while getopts f:o:b: flag
do
    case "${flag}" in
        f) input_f=${OPTARG};;
        o) output_d=${OPTARG};;
        b) database=${OPTARG};;
    esac
done

mkdir -p $output_d

# Search for last common ancestor
function lca(){ cat $@ | sed -e '$!{N;s/^\(.*\).*\n\1.*$/\1\n\1/;D;}' | awk -F ";" '{$NF=""; OFS=";"; print $0}'; return; }

# Perform sequence alignment using usearch
# Search for 500 high identity hits (maxaccepts) in global alignment
usearch -usearch_global $input_f -db $database -id 0.8 -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out $output_d/otus.tax -threads 32 &> $output_d/taxonomy.log

# Extract information from the created 'otus.tax' file and create 'otus.lca'
for i in $(cut -f 1 -d $'\t' $output_d/otus.tax | sort | uniq); do id=$(grep -m 1 -P $i'\t' $output_d/otus.tax | cut -f 3 -d$'\t'); res=$(grep -P $i'\t' $output_d/otus.tax | cut -f 2 -d$'\t' | cut -f 1 -d ' ' --complement | lca); echo -e $i'\t'$res'\t'$id; done > $output_d/otus.lca

