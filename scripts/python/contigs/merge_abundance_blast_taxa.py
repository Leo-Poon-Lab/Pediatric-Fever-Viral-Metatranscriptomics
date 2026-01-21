from sys import argv
import os

sample_name = argv[1]
abundance_file = argv[2]
contigs_file = argv[3]
output_file = argv[4]

# Check if the abundance_file and contigs_file exist and are not empty
if os.path.exists(abundance_file) and os.path.getsize(abundance_file) > 0 and os.path.exists(contigs_file) and os.path.getsize(contigs_file) > 0:

    # Read the abundance data and store it in a dictionary
    abundance_data = {}
    with open(abundance_file, "r") as f:
        header = f.readline().strip()
        for line in f:
            fields = line.strip().split()
            contig = fields[0]
            # Keep only Length, Mapped_Reads, Coverage, and RPM columns (indexes 1, 4, and 6)
            abundance_data[contig] = [fields[1], fields[2], fields[4], fields[6]]

    # Read the contigs data and store it in a dictionary
    contigs_data = {}
    with open(contigs_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            contig = fields[0]
            # Keep only Protein_ID, Identity, Protein_Description, Family, Genus and Species columns
            contigs_data[contig] = [fields[1], fields[2], fields[3], fields[6], fields[8], fields[9], fields[10]]

    # Merge the data and write it to the output file
    with open(output_file, "w") as f:
        # Write the header
        f.write(f"Contig\tLength\tCoverage\tMapped_Reads\tRPM\tProtein_ID\tIdentity\tProtein_Description\tClass\tFamily\tGenus\tSpecies\n")

        # Iterate through the abundance data and merge it with the contigs data
        for contig, abundance_values in abundance_data.items():
            if contig in contigs_data:
                contig_values = contigs_data[contig]
                sample_contig = sample_name + "_" + contig
                f.write(f"{sample_contig}\t{abundance_values[0]}\t{abundance_values[2]}\t{abundance_values[1]}\t{abundance_values[3]}\t{contig_values[0]}\t{contig_values[1]}\t{contig_values[2]}\t{contig_values[3]}\t{contig_values[4]}\t{contig_values[5]}\t{contig_values[6]}\n")
else:
    print(f"Error: {abundance_file} or {contigs_file} does not exist or is empty.")
