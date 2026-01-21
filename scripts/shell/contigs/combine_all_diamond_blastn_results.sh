dir=$1
all_concatenated_merged_blast=$2

files=$(find $dir -type f -name "combined_virus_results.txt")

# Concatenate the first file, including the header, to the output file
head -n 1 $(echo "$files" | head -n 1) > $all_concatenated_merged_blast

# Loop through the rest of the files, remove the header, and append the content to the output file
echo "$files" | tail -n +1 | while read file; do
  # Get the sample name from the file path
  sample_name=$(echo "$file" | awk -F '/' '{print $(NF-2)}')
  # Remove the header, append the sample name to the Contig, and append the content to the output file
  tail -n +2 "$file" | awk -v sample_name="$sample_name" 'BEGIN {FS=OFS="\t"} {gsub(/^/, sample_name "_", $1)} 1' >> $all_concatenated_merged_blast
done