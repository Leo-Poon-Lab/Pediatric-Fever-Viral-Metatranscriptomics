import pandas as pd
from sys import argv
from os import path

work_dir = argv[1]
input_file = path.join(work_dir, 'multiple_taxids.tsv')
output_file = path.join(work_dir, 'split_multiple_taxids.tsv')

# Read Diamond output file into a pandas DataFrame
columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qcovhsp', 'scovhsp', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore','stitle','staxids', 'sskingdoms', 'sscinames']
df = pd.read_csv(input_file, sep='\t', header=None, names=columns)

# Split the staxids column by the semicolon delimiter and create a new DataFrame
staxids_split = df['staxids'].str.split(';').explode().to_frame()

# Join the original DataFrame with the new staxids DataFrame
df_split = df.drop('staxids', axis=1).join(staxids_split)

# Reorder the columns to insert the taxid column in its original order
column_order = columns[:15] + ['staxids'] + columns[16:]
df_split = df_split[column_order]

# Save the new DataFrame with separate rows for each taxid
df_split.to_csv(output_file, sep='\t', index=False, header=None)
