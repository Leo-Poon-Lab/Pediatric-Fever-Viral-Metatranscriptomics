import pandas as pd
import sys
# Read the first file into a DataFrame
input_file = sys.argv[1]
pcr_profile = sys.argv[2]
output_file = sys.argv[3]

df1 = pd.read_csv(input_file, sep='\t')

# Add a new column 'Blastn_Support' based on the numerator of 'Blastn_confirmed'
df1['Blastn_Support'] = df1['Blastn_confirmed'].apply(lambda x: 'N' if int(x.split('/')[0]) == 0 else 'Y')

# Read the Patient_Profile_Specimen_new file into a DataFrame
df2 = pd.read_csv(pcr_profile, sep='\t')

# Merge the two DataFrames on the 'Sample' column
merged_df = pd.merge(df1, df2[['Sample', 'RT-PCR']], on='Sample', how='left')

# Move the RT-PCR column to the second position
cols = list(merged_df.columns)
cols.insert(1, cols.pop(cols.index('RT-PCR')))
merged_df = merged_df[cols]

# Write the merged DataFrame to a new file
merged_df.to_csv(output_file, sep='\t', index=False)
