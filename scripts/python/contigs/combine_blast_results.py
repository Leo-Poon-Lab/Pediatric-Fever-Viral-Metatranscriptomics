import sys

virus_identification_dir = sys.argv[1]
sample_id = sys.argv[2]
blastx_file = sys.argv[3]
blastn_file = sys.argv[4]
combined_file = sys.argv[5]

fasta_file = virus_identification_dir + "/" + sample_id +"/diamond/potential_virus_contigs.fa"

# Read lengths of contigs from the fasta file
contig_lengths = {}
with open(fasta_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            fields = line.strip().split()
            contig_id = fields[0][1:]
            length = int(float(fields[3].split('=')[1]))
            contig_lengths[contig_id] = str(length)

# Read blastx results
blastx_dict = {}
with open(blastx_file, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        contig_id = fields[0]
        sskingdom = fields[16]
        if sskingdom == "Viruses":
            blastx_dict[contig_id] = fields[1:3] + fields[6:8] + fields[14:]
# Read blastn results and select the top match
blastn_dict = {}
with open(blastn_file, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        contig_id = fields[0]
        if contig_id not in blastn_dict:
            blastn_dict[contig_id] = fields[1:5] + [fields[9]] + fields[22:]
        #All the blastn results are stored in a list
        # if contig_id not in blastn_dict:
        #     blastn_dict[contig_id] = []
        # blastn_dict[contig_id].append(fields[1:5] + fields[9:11] + fields[20:22])

# Combine results
combined_results = []
for contig_id in blastx_dict:
    blastx_result = blastx_dict[contig_id]
    blastn_result = blastn_dict.get(contig_id, [])
    contig_length = contig_lengths.get(contig_id, 'NA')
    
    if not blastn_result:
        combined_results.append([contig_id] + [contig_length] + blastx_result)
    else:
        combined_result = [contig_id] + [contig_length] + blastx_result + blastn_result
        combined_results.append(combined_result)

# Headers
header_blastx = ["Contig", "Contig_Length" ,"blastx_SubjectID", "blastx_Identity", "blastx_QueryCoverage", "blastx_SubjectCoverage", "blastx_Title", "blastx_SubjectTaxID", "blastx_SubjectSuperkingdom", "blastx_SubjectSciName"]
header_blastn = [ "blastn_QueryLength", "blastn_SubjectID", "blastn_SubjectLength", "blastn_Identity", "blastn_QueryCoverage", "blastn_SubjectTaxID", "blastn_SubjectTitle"]

# Write combined results with headers
with open(combined_file, "w") as f:
    f.write("\t".join(header_blastx + header_blastn) + "\n")
    for result in combined_results:
        f.write("\t".join(result) + "\n")
