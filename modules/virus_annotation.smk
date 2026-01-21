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
        all_concatenated_file_no_index_hopping_interested = join(abandunce_dir, "quantification_result_interested.txt"),
    output:
        ORF_seqs = join(virus_annotation_dir, "{sample}/ORF_seqs.fa"),
    params:
        min_length = 100,
        contigs_fa = join(virus_annotation_dir, "{sample}/virus_contigs.fa"),

    conda: "base"
    shell:
       """
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
rule select_interesed_proteins:
    input:
        matches = join(virus_annotation_dir, "{sample}/blastp/matches.tsv"),
    output:
        interested_proteins = join(virus_annotation_dir, "{sample}/blastp/interested_proteins.tsv")
    conda: "base"
    params:
        work_dir =("/home/ubuntu/myMetagenomics/Trial/05_virus_annotation/Virus_interested")
    shell:
        """
        #check if the matches are empty
        if [ ! -s {input.matches} ]; then
            touch {output.interested_proteins}
        else
            #从meata data中选择感兴趣的病毒，生成virus_select.csv文件
            (head -1 /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/quantification_result.txt.no.phages.RPM_length_filtered.tsv && grep "Circoviridae" /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/quantification_result.txt.no.phages.RPM_length_filtered.tsv \
            |grep "Arthropods" |sort -k2,2 -nr) > viruse_select.tsv

            # bash virus_select.sh

            #根据病毒的name，到相应的文件夹下找到对应的蛋白质序列
            python select_virus_blastp.py

            #查看修改之后，运行下面的script
            grep -i "replication" viruse_select_blasp.tsv > viruse_select_rep_blasp.tsv

            #查看感兴趣的蛋白质序列，看看是否有RdRp和DNA_pol，如果有，就把这个病毒的蛋白质序列保存下来
            #提取RdRp或DNA pol蛋白质fasta序列
            python extract_and_rename_fasta.py viruse_select_rep_blasp.tsv
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