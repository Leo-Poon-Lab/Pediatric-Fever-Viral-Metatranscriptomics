import csv
from collections import defaultdict
import sys

all_Kraken_Kaiju_combined_blast_results_pcr = sys.argv[1]
all_blast_results_virus_lineages = sys.argv[2]
virus_summary = sys.argv[3]

# Step 1: Read the first file and collect prefixes, samples, RT-PCR status, and Total_Reads
prefix_info = []
with open(all_Kraken_Kaiju_combined_blast_results_pcr, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['Blastn_Support'] == 'Y':
            sample = row['Sample']
            rt_pcr = row['RT-PCR']
            total_reads = int(row['Total_Reads'] )* 2
            taxon_name = row['Taxon_Name']
            # Process taxon name by replacing non-alphanumeric characters with underscores
            processed_taxon = ''.join([c if c.isalnum() else '_' for c in taxon_name])
            prefix = f"{sample}_{processed_taxon}_"
            prefix_info.append((prefix, sample, rt_pcr, total_reads))

# Step 2: Read the blast file and filter rows based on prefixes
filtered_rows = [] 
with open(all_blast_results_virus_lineages, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        query_id = row['Query_ID']
        for prefix, sample, rt_pcr, total_reads in prefix_info:
            if query_id.startswith(prefix):
                # Extract relevant fields
                virus_class = row['Class']
                family = row['Family']
                genus = row['Genus']
                species = row['Species']
                try:
                    perc_identity = float(row['Percentage_Identity'])
                    align_length = int(row['Alignment_Length'])
                except ValueError:
                    # Skip rows with invalid numeric values
                    continue
                filtered_rows.append((sample, rt_pcr, total_reads, virus_class, family, genus, species, perc_identity, align_length))
                break  # Once a prefix matches, move to the next row

# Step 3: Group by sample, RT-PCR, Total_Reads, and taxonomic classification and compute averages
group_stats = defaultdict(lambda: {'count': 0, 'sum_identity': 0.0, 'sum_align': 0})
for entry in filtered_rows:
    sample, rt_pcr, total_reads, v_class, v_family, v_genus, v_species, perc_id, align_len = entry
    key = (sample, rt_pcr, total_reads, v_class, v_family, v_genus, v_species)
    group_stats[key]['count'] += 1
    group_stats[key]['sum_identity'] += perc_id
    group_stats[key]['sum_align'] += align_len

# Prepare the output data
output_data = []
for key in group_stats:
    sample, rt_pcr, total_reads, v_class, v_family, v_genus, v_species = key
    stats = group_stats[key]
    avg_identity = stats['sum_identity'] / stats['count']
    avg_align = stats['sum_align'] / stats['count']
    rpm = stats['count'] / int(total_reads) * 1e6 
    output_data.append([
        sample,
        rt_pcr,
        v_class,
        v_family,
        v_genus,
        v_species,
        total_reads,
        stats['count'],
        rpm,
        f"{avg_identity:.2f}",
        f"{avg_align:.2f}"
    ])

# Step 4: Write the output CSV file
with open(virus_summary, 'w', newline='') as f:
    writer = csv.writer(f,delimiter='\t')
    writer.writerow([
        'Sample', 'RT-PCR', 'Virus_Class', 'Virus_Family', 'Virus_Genus', 'Virus_Species','Total_Reads', 'Number_of_supporting_reads', 'RPM' ,
        'Average_identity', 'Average_alignment_length'
    ])
    writer.writerows(output_data)
