import csv
from sys import argv

quantification_result = argv[1]
ictv_file = argv[2]
# Load ICTV data
ictv_data = {}
with open(ictv_file, 'r') as ictv_file:
    reader = csv.DictReader(ictv_file, delimiter='\t')
    for row in reader:
        family = row['Family']
        species = row['Species']
        genome = row['Genome']
        if family not in ictv_data:
            ictv_data[family] = {}
        ictv_data[family][species] = genome

# Process quantification results and add Genome column
output_rows = []
with open(quantification_result, 'r') as quant_file:

    reader = csv.DictReader(quant_file, delimiter='\t')
    fieldnames = reader.fieldnames + ['Genome']  # Add new column header
    # output_rows.append(dict(zip(fieldnames, fieldnames)))  # Add the updated header to output
    
    for row in reader:
        family = row['Virus_Family']
        species = row['Virus_Species']
        genome = 'Unknown'
        
        if family == 'Totiviridae':
            genome = 'dsRNA'
        elif family in ictv_data:
            if species in ictv_data[family]:
                genome = ictv_data[family][species]
            else:
                # If species not found but family is assigned, use any available genome for that family
                genome = next(iter(ictv_data[family].values()), 'Unknown')
        elif family == 'Unassigned':
            # Search by species if family is unassigned
            for fam, spc_dict in ictv_data.items():
                if species in spc_dict:
                    genome = spc_dict[species]
                    break
                # Additional checks for the Species column if genome is still 'Unknown'
        if genome == 'Unknown':
            species_lower = species.lower()
            if 'totiviridae' in species_lower:
                genome = 'dsRNA'
            elif 'dna' in species_lower:
                genome = 'DNA'
            elif 'rna' in species_lower:
                genome = 'RNA'
            elif 'crucivirus' in species_lower:
                genome = 'ssDNA'
            elif 'cress' in species_lower:
                genome = 'ssDNA'
            elif 'hudisavirus' in species_lower:
                genome = 'DNA'
            elif 'circular' in species_lower:
                genome = 'DNA'
            elif 'reo-like' in species_lower:
                genome = 'dsRNA'
            elif 'uncultured' in species_lower:
                genome = 'dsDNA'
            elif 'pacmanvirus' in species_lower:
                genome = 'DNA'
            elif 'powell' in species_lower:
                genome = 'RNA'
            elif 'chimeric' in species_lower:
                genome = 'DNA'
            elif 'arizlama' in species_lower:
                genome = 'DNA'
            elif 'pinkberry virus' in species_lower:
                genome = 'DNA'
            elif 'pacific flying fox associated' in species_lower:
                genome = 'DNA'
            elif 'biscoe virus' in species_lower:
                genome = 'RNA'
            elif 'bransfield virus' in species_lower:
                genome = 'RNA'
            elif 'dermatophagoides pteronyssinus virus' in species_lower:
                genome = 'RNA'
            elif 'hp38b' in species_lower:
                genome = 'DNA'
            elif 'nioz-uu157' in species_lower:
                genome = 'DNA'
            elif 'virus sp.' in species_lower:
                genome = 'DNA'
            elif 'channeled applesnail virus' in species_lower:
                genome = 'RNA'
            elif 'de lozier virus' in species_lower:
                genome = 'RNA'
            elif 'beihai narna-like virus' in species_lower:
                genome = 'RNA'
            elif 'picorna-like virus' in species_lower:
                genome = 'RNA' 
            elif 'wigfec virus' in species_lower:
                genome = 'DNA'            
            

        row['Genome'] = genome
        output_rows.append(row)


# Write the updated rows back to a new file
# with open('quantification_result_with_genome.txt', 'w', newline='') as outfile:
quantification_result_with_genome_file=argv[3]
with open(quantification_result_with_genome_file, 'w', newline='') as outfile:
    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()  # ✅ This alone adds the header once
    writer.writerows(output_rows)  # Now only contains data rows
