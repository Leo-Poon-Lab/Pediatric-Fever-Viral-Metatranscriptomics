from collections import defaultdict
from sys import argv
import os

def process_kaiju(kaiju_file):
    kaiju_dict = defaultdict(lambda: {'reads': 0, 'taxon_ids': set()})
    with open(kaiju_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            try:
                reads = int(parts[2])
                taxon_id = parts[3]
                taxonomy = parts[4]
            except (IndexError, ValueError):
                continue
            
            # Split taxonomy into components
            # 这里有的taxa是空的，要去掉NA，否则Genus会被当成Species
            taxa = [t.strip() for t in taxonomy.split(';') if t.strip() and t.strip() != 'NA']
            if not taxa or taxa[0] != 'Viruses':
                continue  # Skip if not virus or empty
            
            # Extract ranks explicitly (superkingdom, class, family, genus, species)
            ranks = ['Superkingdom', 'Class', 'Family', 'Genus', 'Species']
            tax_data = {}
            for i, taxon in enumerate(taxa):
                if i < len(ranks):
                    tax_data[ranks[i]] = taxon

            # Use the deepest rank available (species > genus > family > class > superkingdom)
            for rank in reversed(ranks):
                if rank in tax_data:
                    taxon_name = tax_data[rank]
                    key = (taxon_name, rank)
                    kaiju_dict[key]['reads'] += reads
                    kaiju_dict[key]['taxon_ids'].add(taxon_id)
                    break  # Stop at the deepest rank
    return kaiju_dict

def parse_kraken_main_report(main_report_path):
    with open(main_report_path, 'r') as f:
        lines = [line.strip().split() for line in f]
    unclassified = int(lines[0][1])
    classified = int(lines[1][1])
    return unclassified + classified

def process_kraken(kraken_file):
    kraken_dict = defaultdict(lambda: {'reads': 0, 'taxon_ids': set()})
    with open(kraken_file, 'r') as f:
        # next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
            try:
                clade_reads = int(parts[1])
                taxon_id = parts[6]
                rank_code = parts[5]
                genus = parts[9].strip()
                species = parts[10].strip()
            except (IndexError, ValueError):
                continue
            
            # Map Kraken rank codes to our rank names
            rank_map = {
                'S': ('Species', species),
                'G': ('Genus', genus),
                'F': ('Family', parts[8].strip()),
                'C': ('Class', parts[7].strip())
            }

            if rank_code in rank_map:
                rank, taxon_name = rank_map[rank_code]
                if taxon_name != 'NA' and taxon_name:
                    key = (taxon_name, rank)
                    kraken_dict[key]['reads'] += clade_reads
                    kraken_dict[key]['taxon_ids'].add(taxon_id)
    return kraken_dict

def merge_data(kaiju_dict, kraken_dict):
    merged = defaultdict(lambda: {'kaiju_reads': 0, 'kraken_reads': 0, 'taxon_ids': set()})
    # Merge Kaiju entries
    for (name, rank), data in kaiju_dict.items():
        merged[(name, rank)]['kaiju_reads'] += data['reads']
        merged[(name, rank)]['taxon_ids'].update(data['taxon_ids'])
    # Merge Kraken entries
    for (name, rank), data in kraken_dict.items():
        merged[(name, rank)]['kraken_reads'] += data['reads']
        merged[(name, rank)]['taxon_ids'].update(data['taxon_ids'])
    return merged

def write_output(merged_data, sample_name, total_reads, output_file):
    headers = [
        'Sample', 'Taxon_Name', 'Taxon_ID', 'Rank', 
        'Total_Reads', 'Kraken_Reads', 'Kraken_RPM', 'Kaiju_Reads', 
        'Kaiju_RPM'  # New columns
    ]
    
    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for (taxon_name, rank), data in merged_data.items():
            taxon_ids = ','.join(sorted(data['taxon_ids'], key=lambda x: x))
            
            # Calculate RPM values
            kraken_rpm = (data['kraken_reads'] / total_reads) * 1000000 if total_reads != 0 else 0
            kaiju_rpm = (data['kaiju_reads'] / total_reads) * 1000000 if total_reads != 0 else 0
            
            row = [
                sample_name,
                taxon_name,
                taxon_ids,
                rank,
                str(total_reads),
                str(data['kraken_reads']),
                f"{kraken_rpm:.2f}",  # Rounded to 2 decimal places
                str(data['kaiju_reads']),
                f"{kaiju_rpm:.2f}"   # Rounded to 2 decimal places
            ]
            f.write('\t'.join(row) + '\n')

def main():
    # kraken_file = '../../kraken2/PMT9231/PMT9231.kreports.reformatted.virus'
    kraken_file = argv[1]

    # kaiju_file = 'kaiju.out.table'
    kaiju_file = argv[2]

    # main_report_path = '../../kraken2/PMT9231/PMT9231.kreports'
    main_report_path = argv[3]

    sample_name = os.path.basename(main_report_path).split('.')[0]

    total_reads = parse_kraken_main_report(main_report_path)

    kaiju_data = process_kaiju(kaiju_file)
    kraken_data = process_kraken(kraken_file)
    merged = merge_data(kaiju_data, kraken_data)

    output_file = argv[4]
    # output_file = "combined_report.tsv"
    write_output(merged, sample_name,total_reads, output_file)

if __name__ == '__main__':
    main()