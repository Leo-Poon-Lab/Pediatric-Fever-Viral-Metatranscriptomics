import os
import csv
from collections import defaultdict
import re
import sys

def extract_accession(subject_id):
    """Extract accession number from Blastn subject ID"""
    match = re.search(r'gb\|([^|]+)\|', subject_id)
    return match.group(1) if match else subject_id


base_dir = sys.argv[1]
virus_summary_path = sys.argv[2]
combined_blast_path = sys.argv[3]

output_path = "virus_summary_with_contigs.tsv"
output_path = sys.argv[4]

# Initialize structure to store contig data with blast information
contig_data = defaultdict(lambda: {
    'sample': '',
    'length': 0,
    'coverage': 0.0,
    'mapped_reads': 0,
    'rpm': 0.0,
    'protein_id': '',
    'blastx_identities': [],
    'protein_desc': '',
    'class': '',
    'family': '',
    'genus': '',
    'species': '',
    'blastn_subject_ids': set(),
    'blastn_identities': set(),
    'blastn_query_coverages': set(),
    'blastn_titles': set(),
    'blastn_species': set()  # Added for blastn species
})

# Dictionary to map contigs to species
contig_to_species = {}

# First, process the merged_abundance_blast_taxa.txt files
for sample_dir in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample_dir)
    if not os.path.isdir(sample_path):
        continue
    
    contig_file = os.path.join(sample_path, "merged_abundance_blast_taxa.txt")
    if not os.path.exists(contig_file):
        continue
    
    with open(contig_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            contig = row['Contig']
            contig_data[contig]['sample'] = sample_dir
            contig_data[contig]['length'] = int(row['Length'])
            contig_data[contig]['coverage'] = float(row['Coverage'])
            contig_data[contig]['mapped_reads'] = int(row['Mapped_Reads'])
            contig_data[contig]['rpm'] = float(row['RPM'])
            contig_data[contig]['protein_id'] = row['Protein_ID']
            contig_data[contig]['blastx_identities'].append(float(row['Identity']))
            contig_data[contig]['protein_desc'] = row['Protein_Description']
            contig_data[contig]['class'] = row['Class']
            contig_data[contig]['family'] = row['Family']
            contig_data[contig]['genus'] = row['Genus']
            contig_data[contig]['species'] = row['Species']
            
            contig_to_species[contig] = row['Species'].strip()

# Then, add blastn information from the combined file
if os.path.exists(combined_blast_path):
    with open(combined_blast_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            contig = row['Contig']
            if contig in contig_data:
                if 'blastn_SubjectID' in row and row['blastn_SubjectID']:
                    # Extract accession number
                    accession = extract_accession(row['blastn_SubjectID'])
                    contig_data[contig]['blastn_subject_ids'].add(accession)
                    
                    # Add identity if available
                    if 'blastn_Identity' in row and row['blastn_Identity']:
                        contig_data[contig]['blastn_identities'].add(float(row['blastn_Identity']))
                    
                    # Add title if available
                    if 'blastn_SubjectTitle' in row and row['blastn_SubjectTitle']:
                        contig_data[contig]['blastn_titles'].add(row['blastn_SubjectTitle'])
                    
                    # Add query coverage if available
                    if 'blastn_QueryCoverage' in row and row['blastn_QueryCoverage']:
                        contig_data[contig]['blastn_query_coverages'].add(float(row['blastn_QueryCoverage']))
                    
                    # Add blastn species if available
                    if 'blastn_Species' in row and row['blastn_Species']:
                        contig_data[contig]['blastn_species'].add(row['blastn_Species'])

# Now aggregate by sample and species
species_aggregate = defaultdict(lambda: {
    'num_contigs': 0,
    'max_length': 0,
    'max_coverage': 0.0,
    'total_mapped_reads': 0,
    'total_rpm': 0.0,
    'protein_ids': [],
    'blastx_identities': [],
    'protein_descs': [],
    'taxonomy': {},
    'blastn_subject_ids': [],
    'blastn_identities': [],
    'blastn_query_coverages': [],
    'blastn_titles': [],
    'blastn_species': []  # Added for blastn species
})

for contig, data in contig_data.items():
    sample = data['sample']
    species = contig_to_species.get(contig, '')
    if not species:
        continue
    
    key = (sample, species)
    
    species_aggregate[key]['num_contigs'] += 1
    species_aggregate[key]['max_length'] = max(
        species_aggregate[key]['max_length'],
        data['length']
    )
    species_aggregate[key]['max_coverage'] = max(
        species_aggregate[key]['max_coverage'],
        data['coverage']
    )
    species_aggregate[key]['total_mapped_reads'] += data['mapped_reads']
    species_aggregate[key]['total_rpm'] += data['rpm']
    species_aggregate[key]['protein_ids'].append(data['protein_id'])
    species_aggregate[key]['blastx_identities'].extend(data['blastx_identities'])
    species_aggregate[key]['protein_descs'].append(data['protein_desc'])
    
    # Add blastn info if available
    if data['blastn_subject_ids']:
        species_aggregate[key]['blastn_subject_ids'] = list(data['blastn_subject_ids'])
    if data['blastn_identities']:
        species_aggregate[key]['blastn_identities'] = list(data['blastn_identities'])
    if data['blastn_titles']:
        species_aggregate[key]['blastn_titles'] = list(data['blastn_titles'])
    if data['blastn_query_coverages']:
        species_aggregate[key]['blastn_query_coverages'] = list(data['blastn_query_coverages'])
    if data['blastn_species']:  # Added for blastn species
        species_aggregate[key]['blastn_species'] = list(data['blastn_species'])
    
    # Store taxonomy from first occurrence
    if species_aggregate[key]['num_contigs'] == 1:
        species_aggregate[key]['taxonomy'] = {
            'Class': data['class'],
            'Family': data['family'],
            'Genus': data['genus']
        }

# Dictionary to store PCR values from virus_summary.tsv
pcr_dict = {}

# Merge with virus_summary
with open(virus_summary_path, 'r') as f_in, open(output_path, 'w') as f_out:
    reader = csv.DictReader(f_in, delimiter='\t')
    fieldnames = reader.fieldnames + [
        'Number_of_contigs',
        'Length',
        'Coverage',
        'Mapped_Reads',
        'RPM_mapping_contigs',
        'Protein_ID',
        'Blastx_Identity',
        'Protein_Description',
        'Blastn_SubjectID',
        'Blastn_Identity',
        'Blastn_QueryCoverage',
        'Blastn_Title',
        'Blastn_Species'  # Added for blastn species
    ]
    writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    
    processed_keys = set()
    
    # Process existing records
    for row in reader:
        key = (row['Sample'], row['Virus_Species'])
        processed_keys.add(key)
        
        # Store the PCR value in the dictionary
        pcr_dict[key[0]] = row['PCR']

        if key in species_aggregate:
            data = species_aggregate[key]
            row.update({
                'Number_of_contigs': data['num_contigs'],
                'Length': data['max_length'],
                'Coverage': data['max_coverage'],
                'Mapped_Reads': data['total_mapped_reads'],
                'RPM_mapping_contigs': data['total_rpm'],
                'Protein_ID': ';'.join(data['protein_ids']),
                'Blastx_Identity': ';'.join(map(str, data['blastx_identities'])),
                'Protein_Description': ';'.join(data['protein_descs']),
                'Blastn_SubjectID': ';'.join(data['blastn_subject_ids']),
                'Blastn_Identity': ';'.join(map(str, data['blastn_identities'])),
                'Blastn_QueryCoverage': ';'.join(map(str, data['blastn_query_coverages'])),
                'Blastn_Title': ';'.join(data['blastn_titles']),
                'Blastn_Species': ';'.join(data['blastn_species'])  # Added for blastn species
            })
        else:
            # Add empty values if no contig data
            row.update({
                'Number_of_contigs': '',
                'Length': '',
                'Coverage': '',
                'Mapped_Reads': '',
                'RPM_mapping_contigs': '',
                'Protein_ID': '',
                'Blastx_Identity': '',
                'Protein_Description': '',
                'Blastn_SubjectID': '',
                'Blastn_Identity': '',
                'Blastn_QueryCoverage': '',
                'Blastn_Title': '',
                'Blastn_Species': ''  # Added for blastn species
            })
        
        writer.writerow(row)

    # Add new records from contig data that weren't in virus_summary
    for key in species_aggregate:
        if key not in processed_keys:
            data = species_aggregate[key]
            taxonomy = data['taxonomy']

            new_row = {
                'Sample': key[0],
                'PCR': pcr_dict.get(key[0], ''),
                'Virus_Class': taxonomy.get('Class', ''),
                'Virus_Family': taxonomy.get('Family', ''),
                'Virus_Genus': taxonomy.get('Genus', ''),
                'Virus_Species': key[1],
                'Total_Reads': '',
                'Number_of_supporting_reads': '',
                'RPM': '',
                'Average_identity': '',
                'Average_alignment_length': '',
                'Number_of_contigs': data['num_contigs'],
                'Length': data['max_length'],
                'Coverage': data['max_coverage'],
                'Mapped_Reads': data['total_mapped_reads'],
                'RPM_mapping_contigs': data['total_rpm'],
                'Protein_ID': ';'.join(data['protein_ids']),
                'Blastx_Identity': data['blastx_identities'],
                'Protein_Description': ';'.join(data['protein_descs']),
                'Blastn_SubjectID': ';'.join(data['blastn_subject_ids']),
                'Blastn_Identity': ';'.join(map(str, data['blastn_identities'])),
                'Blastn_QueryCoverage': ';'.join(map(str, data['blastn_query_coverages'])),
                'Blastn_Title': ';'.join(data['blastn_titles']),
                'Blastn_Species': ';'.join(data['blastn_species'])  # Added for blastn species
            }
            writer.writerow(new_row)

print(f"Contig and blastn information merged successfully! Output written to {output_path}")