rule filter_phage_rpm:
    input:
        all_concatenated_file = join(abandunce_dir, "quantification_result.txt"),
    output:
        all_concatenated_file_filtered_phages= join(abandunce_dir, "quantification_result.txt.no.phages"),
        all_concatenated_file_filtered_phages_rpm= join(abandunce_dir, "quantification_result.txt.no.phages.RPM"),
    conda: "base"
    shell:
        """
        # Filter phages
        grep -vE "Bacteria|phage" {input.all_concatenated_file} > {output.all_concatenated_file_filtered_phages}
        
        # Filter RPM < 1
        awk -F"," '$5 >= 1' {input.all_concatenated_file_no_index_hopping}.no.phages > {input.all_concatenated_file_filtered_phages_rpm}
        """

rule orfFinder:
    input:
        # all_concatenated_file_no_index_hopping = join(abandunce_dir, "quantification_result.txt"),
        all_concatenated_file_no_index_hopping_interested = join(abandunce_dir, "quantification_result_interested.txt"),
    output:
        ORF_seqs = join(virus_annotation_dir, "{sample}/ORF_seqs.fa"),
    params:
        min_length = 100,
        contigs_fa = join(virus_annotation_dir, "{sample}/virus_contigs.fa"),

    conda: "base"
    shell:
       """
        # # Filter phages
        # grep -vE "Bacteria|phage" {input.all_concatenated_file_no_index_hopping} > {input.all_concatenated_file_no_index_hopping}.no.phages

        # # Filter RPM < 1
        # awk -F"," '$5 >= 1' {input.all_concatenated_file_no_index_hopping}.no.phages > {input.all_concatenated_file_no_index_hopping}.no.phages.RPM

        # #generate fasta file
        # rm {params.contigs_fa}
        
        # python /home/ubuntu/Tools/Scripts/generate_contigs_file_using_results.py {input.all_concatenated_file_no_index_hopping}.no.phages.RPM
        python /home/ubuntu/Tools/Scripts/generate_contigs_file_using_results.py {input.all_concatenated_file_no_index_hopping_interested}

        if [ ! -s {params.contigs_fa} ]; then
            touch {output.ORF_seqs}
        else
            /home/ubuntu/Tools/ORFfinder/ORFfinder -in {params.contigs_fa} -ml {params.min_length} -out {output.ORF_seqs}
        fi
        """

checkpoint blastp:
    input:
        ORF_seqs = join(virus_annotation_dir, "{sample}/ORF_seqs.fa"),
    output:
        matches = join(virus_annotation_dir, "{sample}/blastp/matches.tsv")
    threads: 10
    params:
        database = "/home/ubuntu/Tools/Diamond/nr_tax_full.dmnd",
        evalue = 1e-5
    benchmark:
        join(virus_annotation_dir, "{sample}/blastp/{sample}.blastp.benchmark.txt")
    conda: "diamond"
    shell:
        """
        #check if the contigs are empty
        if [ ! -s {input.ORF_seqs} ]; then
            touch {output.matches}
        else
            diamond blastp -q {input.ORF_seqs} -d {params.database} -o {output.matches} \
            --very-sensitive --evalue {params.evalue} --max-target-seqs 1 --threads {threads} \
            -f 6 qseqid sseqid pident length mismatch gapopen qcovhsp scovhsp qstart qend sstart send evalue bitscore stitle staxids sskingdoms sscinames
        fi
        """

# checkpoint diamond:
#     input:
#         contigs = join(assembly_dir, "{sample}/megahit/final.contigs.fa"),
#         # contigs = join(assembly_dir, "{sample}/filtered_contigs/filtered.final.small.contigs.fa"),
#     output:
#         matches = join(virus_identification_dir, "{sample}/diamond/matches.tsv")
#     threads: 10
#     params:
#         database = "/home/ubuntu/Tools/Diamond/nr_tax_full.dmnd",
#         evalue = 1e-5
#     benchmark:
#         join(virus_identification_dir, "{sample}/diamond/{sample}.diamond.benchmark.txt")
#     conda: "diamond"
#     shell:
#         """
#         #check if the contigs are empty
#         if [ ! -s {input.contigs} ]; then
#             touch {output.matches}
#         else
#             diamond blastx -q {input.contigs} -d {params.database} -o {output.matches} \
#             --more-sensitive --evalue {params.evalue} --max-target-seqs 1 --threads {threads} \
#             -f 6 qseqid sseqid pident length mismatch gapopen qcovhsp scovhsp qstart qend sstart send evalue bitscore stitle staxids sskingdoms sscinames
#         fi
#         """

# # select kingdom virus as potential virus
# rule identify_potential_virus:
#     input:
#         matches = join(virus_identification_dir, "{sample}/diamond/matches.tsv"),
#         contigs = join(assembly_dir, "{sample}/megahit/final.contigs.fa"),
#     output:
#         virus_matches = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs_all.tsv"),
#         virus_matches_contigs = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs.fa")
#     params:
#         matches_dir = join(virus_identification_dir, "{sample}/diamond")
#     conda: "base"
#     shell:
#         """
#         #check if the matches are empty
#         if [ ! -s {input.matches} ]; then
#             touch {output.virus_matches}
#             touch {output.virus_matches_contigs}
#         else
#             # generate formal taxonomic table
#             bash /home/ubuntu/Tools/Scripts/taxonkit_retrive_taxa_from_taxaid.sh {params.matches_dir} 
#             if grep -qi "Virus" {output.virus_matches}; then
#                 seqkit grep -f <(cut -f1 {output.virus_matches}) {input.contigs} > {output.virus_matches_contigs}
#             else
#                 touch {output.virus_matches}
#             fi
#         fi
        
#         """

# # rule identify_unhit_contigs:
# #     input:
# #         contigs = join(assembly_dir, "{sample}/filtered_contigs/filtered.final.contigs.fa"),
# #         matches = join(virus_identification_dir, "{sample}/diamond/matches.tsv")
# #     output:
# #         unhit_contigs = join(virus_identification_dir, "{sample}/diamond/missing_contigs.fa")
# #     params:
# #         contig_ids = join (virus_identification_dir, "{sample}/diamond/contig_ids.txt"),
# #         matched_ids = join (virus_identification_dir, "{sample}/diamond/matched_ids.txt"),
# #         missing_ids = join (virus_identification_dir, "{sample}/diamond/missing_ids.txt"),
# #         matches_dir = join(virus_identification_dir, "{sample}/diamond")
# #     conda: "base"
# #     shell:
# #         """
# #         grep ">" {input.contigs}|  sed 's/^>//g' | cut -d' ' -f1 |sort|uniq > {params.contig_ids}
# #         cut -f1 {input.matches}| sort | uniq  > {params.matched_ids}
# #         comm -23 {params.contig_ids} {params.matched_ids} > {params.missing_ids}
# #         seqkit grep -f {params.missing_ids} {input.contigs} > {output.unhit_contigs}
#         # """



# rule anlysis_blastn_blastx_virus:
#     input:
#         virus_matches_diamond = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs_all.tsv"),
#         virus_matches_blastn = join(virus_identification_dir, "{sample}/blastn/matches.tsv"),
#     output:
#         no_conflict_contigs = join(virus_identification_dir, "{sample}/comprenhensive_anlysis/no_conflict_contigs.txt")
#     params:
#         work_dir = join(virus_identification_dir, "{sample}/comprenhensive_anlysis"),
#         blastn_dir = join(virus_identification_dir, "{sample}/blastn"),
#     conda: "base"
#     shell:
#         """
#         bash /home/ubuntu/Tools/Scripts/anlysis_blastn_blastx.sh {params.work_dir} {params.blastn_dir} {input.virus_matches_blastn} \
#         {input.virus_matches_diamond} {output.no_conflict_contigs}
        
#         """

# rule diamond_rdrp_dna_pol:
#     input:
#         virus_matches_contigs = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs.fa")
#         # contigs = join(assembly_dir, "{sample}/filtered_contigs/filtered.final.small.contigs.fa"),
#     output:
#         matches = join(virus_identification_dir, "{sample}/diamond_rdrp_dna/matches_rdrp_dna.tsv")
#     threads: 20
#     params:
#         database = "/home/ubuntu/Tools/Diamond/virus/RdRp_DNA_pol/RdRp_DNA_pol_db.dmnd",
#         evalue = 1e-5
#     benchmark:
#         join(virus_identification_dir, "{sample}/diamond_rdrp_dna/{sample}.diamond.benchmark.txt")
#     conda: "diamond"
#     shell:
#         """
#         #check if the contigs are empty
#         if [ ! -s {input.virus_matches_contigs} ]; then
#             touch {output.matches}
#         else
#             diamond blastx -q {input.virus_matches_contigs} -d {params.database} -o {output.matches} \
#             --very-sensitive --evalue {params.evalue} --max-target-seqs 1 --threads {threads} \
#             -f 6 qseqid sseqid pident length mismatch gapopen qcovhsp scovhsp qstart qend sstart send evalue bitscore stitle staxids sskingdoms sscinames
#         fi
#         """