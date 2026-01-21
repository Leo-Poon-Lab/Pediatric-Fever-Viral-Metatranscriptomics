#merge_virus_abundance_taxa_host
# Get the list of files to concatenate
dir=$1
all_concatenated_merged_abundance_blast_taxa=$2

files=$(find $dir -type f -name "merged_abundance_blast_taxa_grouped_by_species.txt")

# Concatenate the first file, including the header, to the output file
cat $(echo "$files" | head -n 1) > $all_concatenated_merged_abundance_blast_taxa

# Loop through the rest of the files, remove the header, and append the content to the output file
echo "$files" | tail -n +2 | while read file; do
  tail -n +2 "$file" >> $all_concatenated_merged_abundance_blast_taxa
done

#cut -f9 all_concatenated_merged_abundance_blast_taxa.txt > species.txt
(head -n 1 $all_concatenated_merged_abundance_blast_taxa|cut -f9,10,11,12; cut -f9,10,11,12 $all_concatenated_merged_abundance_blast_taxa|tail -n +2 | sort -u) > $dir/class_family_genus_species.txt
