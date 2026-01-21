import csv
from sys import argv
# This script takes the output of the merge_abundance_blast_taxa.py script and groups the rows by species, sums the RPM values, and keeps the longest length for each species.
input_file = argv[1]
output_file = argv[2]

def read_input_file(input_file):
    with open(input_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        data = [row for row in reader]
    return headers, data

def group_by_species(data):
    species_index = -1
    for i, header in enumerate(headers):
        if header == "Species":
            species_index = i
            break

    species_groups = {}
    for row in data:
        species = row[species_index]
        if species not in species_groups:
            species_groups[species] = []
        species_groups[species].append(row)

    return species_groups

def sum_rpm_and_longest_length(species_groups, headers):
    length_index = headers.index("Length")
    rpm_index = headers.index("RPM")
    mapped_reads_index = headers.index("Mapped_Reads")

    result = []
    for species, group in species_groups.items():
        filtered_group = [row for row in group if float(row[rpm_index]) > 0]
        # filtered_group = [row for row in group if float(row[rpm_index]) > 1]
        if filtered_group:
            max_length = max(int(row[length_index]) for row in filtered_group)
            total_rpm = sum(float(row[rpm_index]) for row in filtered_group)
            total_mapped_reads = sum(int(row[mapped_reads_index]) for row in filtered_group)

            representative_row = [row for row in filtered_group if int(row[length_index]) == max_length][0]
            representative_row[rpm_index] = total_rpm
            representative_row[mapped_reads_index] = total_mapped_reads
            result.append(representative_row)

    return result

def write_output_file(headers, data, output_file):
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        for row in data:
            writer.writerow(row)

headers, data = read_input_file(input_file)
species_groups = group_by_species(data)
result = sum_rpm_and_longest_length(species_groups, headers)
write_output_file(headers, result, output_file)
