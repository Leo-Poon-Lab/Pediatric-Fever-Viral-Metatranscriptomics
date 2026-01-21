#integration of all blast files of samples
reads_taxa_kraken_kaiju_dir=$1
all_blast_results_virus_lineages=$2

find $reads_taxa_kraken_kaiju_dir -name *_blast_with_lineage.tsv | xargs awk 'FNR==1 && NR!=1 {next;}{print}'> $reads_taxa_kraken_kaiju_dir/all_blast_results.tsv
awk -F'\t' 'NR == 1 || tolower($26) ~ /virus/' $reads_taxa_kraken_kaiju_dir/all_blast_results.tsv > $reads_taxa_kraken_kaiju_dir/all_blast_results_virus.tsv
cut -f23 $reads_taxa_kraken_kaiju_dir/all_blast_results_virus.tsv | sed '1d' | taxonkit reformat -I 1 -r Unassigned -f "{c}\t{f}\t{g}\t{s}" \
            | sed '1i\TaxID\tClass\tFamily\tGenus\tSpecies' \
            | paste $reads_taxa_kraken_kaiju_dir/all_blast_results_virus.tsv - > $all_blast_results_virus_lineages

rm  $reads_taxa_kraken_kaiju_dir/all_blast_results.tsv
rm $reads_taxa_kraken_kaiju_dir/all_blast_results_virus.tsv
