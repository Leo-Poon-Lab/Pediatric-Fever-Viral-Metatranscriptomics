import csv
from sys import argv

def reformat_krakenuniq(input_path, output_path):
    # Define the hierarchy of taxonomic ranks in order
    rank_order = ['D', 'K', 'P', 'C', 'O', 'F', 'G', 'S']
    current_path = {}  # Tracks the current path of taxonomic ranks
    entries = []       # Collects relevant entries for output
    
    with open(input_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 8:
                continue
            
            rank_code = parts[5]
            taxon_name = ' '.join(parts[7:])
            rank_char = rank_code[0] if rank_code else ''
            
            # Update current_path for recognized ranks and clear lower ranks
            if rank_char in rank_order:
                try:
                    idx = rank_order.index(rank_char)
                except ValueError:
                    continue  # Shouldn't happen due to previous check
                
                # Reset all lower ranks
                for lower_rank in rank_order[idx+1:]:
                    if lower_rank in current_path:
                        del current_path[lower_rank]
                # Update current rank
                current_path[rank_char] = taxon_name
            
            # Process entries of interest (C, F, G, S) under Viruses
            if rank_char in ['C', 'F', 'G', 'S']:
                # Check if current domain (D) contains 'viruses' (case-insensitive)
                current_class = current_path.get('C', '').lower()
                current_family = current_path.get('F', '').lower()
                current_genus = current_path.get('G', '').lower()
                current_species = taxon_name.lower() if rank_char == 'S' else ''

                if any('virus' in term for term in [current_class, current_family, current_genus, current_species]):
                    # Prepare entry with current taxonomy path
                    entry = {
                        'fields': parts[:7],
                        'rank': rank_code,
                        'class': current_path.get('C', 'NA'),
                        'family': current_path.get('F', 'NA'),
                        'genus': current_path.get('G', 'NA'),
                        'species': taxon_name if rank_char == 'S' else 'NA'
                    }
                    entries.append(entry)

    # Filter entries to exclude parent taxa when child exists
    filtered = []
    seen_parents = set()
    # Iterate in reverse to find child entries first
    for entry in reversed(entries):
        if entry['rank'] == 'S':
            filtered.append(entry)
            # 标记父级Genus需要排除
            seen_parents.add(entry['genus'])
        elif entry['rank'] == "G" and entry['genus'] not in seen_parents:
            filtered.append(entry)


    # Write the filtered results to the output file
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(['Percentage', 'CladeReads', 'TaxonReads', 'minimizers', 
                         'distinct_minimizers', 'Rank', 'TaxID', 'Class', 'Family', 
                         'Genus', 'Species'])
        # Output in original order, filtered to remove redundant parents
        for entry in reversed(filtered):
            row = entry['fields'] + [
                entry['class'],
                entry['family'],
                entry['genus'],
                entry['species']
            ]
            writer.writerow(row)


# Example usage
kreports = argv[1]
reformatted_virus_kreports = argv[2]

reformat_krakenuniq(kreports, reformatted_virus_kreports)