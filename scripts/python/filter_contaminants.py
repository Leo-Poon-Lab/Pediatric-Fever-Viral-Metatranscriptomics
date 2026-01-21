#!/usr/bin/env python3
from sys import argv
def filter_virus_summary(virus_file, contamination_file, retained_file):
    # Load contamination descriptions
    contaminants = set()
    with open(contamination_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 1:
                contaminants.add(parts[0].strip())

    # Process virus summary
    with open(virus_file, 'r') as fin, \
        open(retained_file, 'w') as f_retained:

        # Write headers
        header = next(fin)
        f_retained.write(header)

        # Process each entry
        for line in fin:
            cols = line.strip().split('\t')
            # Check relevant fields for contamination
            virus_fields = {
                cols[2].strip(),  # Virus_Class
                cols[3].strip(),  # Virus_Family
                cols[4].strip(),  # Virus_Genus
                cols[5].strip(),   # Virus_Species
                cols[24].strip()   # Blastn_Title
            }

            if virus_fields & contaminants:  # Set intersection
                continue
            else:
                f_retained.write(line)

if __name__ == "__main__":

    input_file = argv[1]
    contamination_file = argv[2]
    retained_file = argv[3]
    
    filter_virus_summary(input_file, contamination_file, retained_file)
