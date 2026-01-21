import pandas as pd
from sys import argv
# # Read the CSV files
quantification_file = argv[1]
metagenomic_samples_file = argv[2]
filtered_output = argv[3]
removed_output = argv[4]  # New output file for filtered rows

# Read data
quantification_result = pd.read_csv(quantification_file, sep="\t")
metagenomic_samples = pd.read_csv(metagenomic_samples_file, sep= "\t")

# Create sample-group mapping
sample_group_dict = dict(zip(metagenomic_samples['Sample'], metagenomic_samples['Group']))

# Add Group column
def extract_group(contig):
    parts = contig.split('_')
    if parts[0] == "PicoInput":
        return "_".join(parts[:2])
    return parts[0]

quantification_result['Group'] = quantification_result['Contig'].apply(extract_group).map(sample_group_dict)

# Calculate maximum reads per group-species combination
quantification_result['Max_in_Group_Species'] = quantification_result.groupby(
    ['Group', 'Species']
)['Mapped_Reads'].transform('max')

# Create filter mask
filter_mask = quantification_result['Mapped_Reads'] >= 0.001 * quantification_result['Max_in_Group_Species']

# Split data
filtered_result = quantification_result[filter_mask].drop(columns=['Max_in_Group_Species'])
removed_rows = quantification_result[~filter_mask].drop(columns=['Max_in_Group_Species'])

# Save results
filtered_result.to_csv(filtered_output, index=False, sep="\t")
removed_rows.to_csv(removed_output, index=False, sep="\t")
