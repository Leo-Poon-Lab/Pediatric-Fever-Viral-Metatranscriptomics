all_combined_virus_results=$1
all_combined_virus_results_lineages=$2
all_combined_virus_results_no_phages=$3
all_combined_virus_results_no_phages_lineages=$4
cut -f16 $all_combined_virus_results | sed '1d' | sed 's/^$/Unassigned/' | taxonkit reformat -I 1 -r Unassigned -f "{c}\t{f}\t{g}\t{s}" |  sed '1i\blastn_TaxID\tblastn_Class\tblastn_Family\tblastn_Genus\tblastn_Species' | paste $all_combined_virus_results - > $all_combined_virus_results_lineages

cut -f15 $all_combined_virus_results_no_phages | sed '1d' | sed 's/^$/Unassigned/' | taxonkit reformat -I 1 -r Unassigned -f "{c}\t{f}\t{g}\t{s}" |  sed '1i\blastn_TaxID\tblastn_Class\tblastn_Family\tblastn_Genus\tblastn_Species' | paste $all_combined_virus_results_no_phages - > $all_combined_virus_results_no_phages_lineages
