from sys import argv
from os import path

work_dir = argv[1]
input_file = path.join(work_dir, 'multiple_taxids_taxa.txt')
output_file = path.join(work_dir, 'potential_virus_contigs_2.tsv')
# Read the file
with open(input_file, 'r') as file:
    lines = file.readlines()

# Dictionary to store the sequences
sequences = {}

# Iterate over each line in the file
for line in lines:
    # Split the line into columns
    columns = line.strip().split('\t')

    # Extract the sequence ID and domain
    sequence_id = columns[0]
    domain = columns[4]
    virus_class = columns[6]
    # Check if the row belongs to Virus
    if domain == 'Viruses':
        # Check if the sequence ID is already in the dictionary
        if sequence_id in sequences:
            # Check if the current row has a different class than the existing row
            if virus_class != 'Caudoviricetes':
                # Update the row
                sequences[sequence_id] = columns
        else:
            # Add the row to the dictionary
            sequences[sequence_id] = columns 

# Write the filtered sequences to a new file
with open(output_file, 'w') as file:
    file.write('\n'.join(['\t'.join(columns) for columns in sequences.values()]))
