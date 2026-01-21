#!/bin/bash

# Input and output file paths
# input_file="kraken_kaiju_combined_report.tsv"
# output_file="kraken_kaiju_combined_report_with_class.tsv"
# temp_file="temp_taxids.txt"
# temp_taxonkit_output="temp_taxonkit_output.txt"

input_file=$1
output_file=$2
temp_file=$3
temp_taxonkit_output=$4
unique_taxids=$5

# Extract tax IDs and process them
tail -n +2 "$input_file" | cut -f3 > "$temp_file"

# Create the output file and add the header
head -n 1 "$input_file" | awk '{print $0"\tClass"}' > "$output_file"

# Function to process tax IDs
process_taxids() {
    local taxids="$1"
    IFS=',' read -r -a taxid_array <<< "$taxids"

    local classes=()
    for taxid in "${taxid_array[@]}"; do
        class=$(grep "^$taxid" "$temp_taxonkit_output" | cut -f2)
        classes+=("$class")
    done

    unique_classes=($(printf "%s\n" "${classes[@]}" | sort -u))
    if [ "${#unique_classes[@]}" -eq 1 ]; then
        echo "${unique_classes[0]}"
    else
        echo "$(IFS=,; echo "${unique_classes[*]}")"
    fi
}

# Create a file with all unique tax IDs
tr ',' '\n' < "$temp_file" | sort -u > $unique_taxids

# Run taxonkit on the file with all unique tax IDs
# 形成一个字典
# 10379   Herviviricetes
# 11855   Revtraviricetes
# 12058   Pisoniviricetes
taxonkit reformat -I 1 -r Unassigned -f "{c}" < $unique_taxids > "$temp_taxonkit_output"

# Process each line in the original temp file and write to output file
while IFS= read -r line; do
    taxids=$(echo "$line" | cut -f3)
    tax_class=$(process_taxids "$taxids")
    echo -e "$line\t$tax_class" >> "$output_file"
done < <(tail -n +2 "$input_file")

# Clean up
rm "$temp_file"
rm "$unique_taxids"
# rm "$temp_taxonkit_output"