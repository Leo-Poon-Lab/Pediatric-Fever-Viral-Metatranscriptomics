#!/bin/bash

# Directory containing BLAST files
BLAST_DIR=$1

# Loop through all BLAST files in the directory
for file in "$BLAST_DIR"/*.combined_blast.tsv; do
    # Define the output file name
    output_file="${file%.tsv}_with_lineage.tsv"
    # Extract Subject_Taxonomy_ID (23rd column), get lineage, and paste back into the file
    cut -f23 "$file" | \
    taxonkit lineage --show-lineage-taxids|sed 1d|sed '1iSubject_Taxonomy_ID\tLineage\tTaxIDs' | \
    paste "$file" - > "$output_file"

    # echo "Done processing $file"
done
