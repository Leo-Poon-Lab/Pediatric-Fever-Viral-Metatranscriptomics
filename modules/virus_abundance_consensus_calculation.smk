rule coverage_reads_reference_alignment_pipeline:
    input:
        rmhost_r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
        rmhost_r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz"),
        virus_contigs = join(abandunce_consensus_dir, "{sample}/samples_reference_info.tsv")
    output:
        aligned_reads = join(abandunce_consensus_dir, "{sample}/virus_aligned.bam"),
    conda: "base"
    threads: 5
    params:
        index_dir = join(abandunce_consensus_dir, "{sample}/viral_contigs_index"),
        index = join(abandunce_consensus_dir, "{sample}/viral_contigs_index/viral_contigs_index"),
        sensitivity = "--very-sensitive"
    shell:
        """
            bash /home/ubuntu/myMetagenomics/Trial/04_abundance_consensus/process_samples_pipeline.sh <Sample> <Species>
            # Rscript /home/ubuntu/myMetagenomics/Trial/04_abundance_consensus/plot_coverage.R Shenzhen_dicistro-like_virus.MT757492.1_coverage.txt Shenzhen_dicistro-like_virus.MT757492.1.gff Shenzhen_dicistro-like_virus.MT757492.1_coverage_plot.pdf
        """

rule ivar_consensus:
    input:
        rmhost_r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
        rmhost_r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz"),
        virus_contigs = join(abandunce_consensus_dir, "{sample}/samples_reference_info.tsv")
    output:
        aligned_reads = join(abandunce_consensus_dir, "{sample}/virus_aligned.bam"),
    conda: "phylogeny"
    threads: 5
    params:
        index_dir = join(abandunce_consensus_dir, "{sample}/viral_contigs_index"),
        index = join(abandunce_consensus_dir, "{sample}/viral_contigs_index/viral_contigs_index"),
        sensitivity = "--very-sensitive"
    shell:
        """
            samtools mpileup -A -d 0 -Q 0 Rhinovirus_B.MN306025.1_aligned.bam | ivar consensus -p Respirovirus_pneumoniae_consensus -t 0.75 -m 3
            echo "N% = $(grep -v '^>' Respirovirus_pneumoniae_consensus.fa | tr -d '\n' | awk '{n=gsub(/N|n/, ""); print n*100/length}')%"
            echo "N = $(grep -v '^>' Respirovirus_pneumoniae_consensus.fa | tr -d '\n' | awk '{n=gsub(/N|n/, ""); print n}')"

        """

# rule coverage_stat:
#     input:
#         aligned_reads = join(abandunce_dir, "{sample}/virus_aligned.bam")
#     output:
#         coverage_log = join(abandunce_dir, "{sample}/weesam_virus_coverage.txt"),
#         coverage_sorted = join(abandunce_dir, "{sample}/coverage_sorted.tsv"),
#     conda: "base"
#     threads: 10
#     params:
#         file_dir = join(abandunce_dir, "{sample}"),
#         wee_sam_html = join(abandunce_dir, "{sample}/wee_sam.html"),
#     shell:
#         """
#         if [ ! -s {input.aligned_reads} ]; then
#             touch {output.coverage_log}
#             touch {output.coverage_sorted}
#             exit 0
#         else        
#             /home/ubuntu/Tools/weeSAM/weeSAM-master/weeSAM --bam {input.aligned_reads} \
#             --out {output.coverage_log} --html {params.wee_sam_html}  --overwrite  --pretty > {params.file_dir}/weesam.log 2>&1

#             #截取有用列，挑选符合阈值的记录
#             #默认coverage>=50% & depth>=1，以depth为主要排序，其次是coverage
#             #暂时不进行过滤
#             cat {output.coverage_log} | cut -f1,2,3,4,5,8 | awk -v FS="\t" -v OFS="\t" '{{if(($5>=50)&&($6>=1)) print $0}}' | sort -t $'\t' -k6nr -k5nr > {params.file_dir}/coverage_depth_sorted.tsv
#             cat {output.coverage_log} | cut -f1,2,3,4,5,8 | awk -v FS="\t" -v OFS="\t" '{{if($5>=0) print $0}}' | sort -t $'\t' -k5nr -k6nr > {params.file_dir}/coverage_sorted.tsv
#         fi

#         """

# rule virus_abundance_estimation:
#     input:
#         aligned_reads = join(abandunce_dir, "{sample}/virus_aligned.bam"),
#         virus_coverage_file = join(abandunce_dir, "{sample}/coverage_sorted.tsv"),
#     output:
#         virus_abundance = join(abandunce_dir, "{sample}/virus_abundance.txt"),
#     conda: "base"
#     params:
#         file_dir = join(abandunce_dir, "{sample}"),
#         sample_name = "{sample}"
#     shell:
#         """
#         total_reads=$(samtools view -c {input.aligned_reads})
#         echo -e "Contig\tLength\tMapped_Reads\tBreadth\tCoverage\tAvg_Depth\tRPM" > {output.virus_abundance}

#         awk -v total_reads=$total_reads \
#         'BEGIN {{ OFS="\t" }} {{ rpm = ($3 / total_reads) * 1000000; \
#         printf("%s\t%s\t%s\t%s\t%s\t%s\t%.2f\\n", $1, $2, $3, $4, $5, $6, rpm) }}' \
#         {input.virus_coverage_file} >> {output.virus_abundance}
        
#         """

# rule merge_virus_abundance_taxa:
#     input:
#         virus_abundance = join(abandunce_dir, "{sample}/virus_abundance.txt"),
#         virus_matches = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs_all_no_phage.tsv"),
#     output:
#         virus_abundance_taxa = join(abandunce_dir, "{sample}/merged_abundance_blast_taxa.txt"),
#         virus_abundance_taxa_grouped = join(abandunce_dir, "{sample}/merged_abundance_blast_taxa_grouped_by_species.txt"),
#     conda: "base"
#     params:
#         sample_name = "{sample}"
#     shell:
#         """ 
#         if [ ! -s {input.virus_matches} ]; then
#             touch {output.virus_abundance_taxa}
#             touch {output.virus_abundance_taxa_grouped}
#             exit 0
#         else  
# 		    #生成第一个contigs的注释文件merged_abundance_blast_taxa.txt，列名有：Contig  Length  Coverage        Mapped_Reads    RPM     Protein_ID      Identity        Protein_Description     Taxonomy Info...
#             python /home/ubuntu/Tools/Scripts/merge_abundance_blast_taxa.py {params.sample_name} {input.virus_abundance} {input.virus_matches} {output.virus_abundance_taxa}
            
#             #这是将同一个样本中同物种species的contigs只用一个最长的contig来进行表示，并且将它们对应的RPM进行相加
#             python /home/ubuntu/Tools/Scripts/sum_rpm_and_group_species.py {output.virus_abundance_taxa} {output.virus_abundance_taxa_grouped}

#             if [ ! -s "/home/ubuntu/myMetagenomics/Trial/04_abundance/all_concatenated_merged_abundance_blast_taxa.txt" ]; then
#                 touch /home/ubuntu/myMetagenomics/Trial/04_abundance/all_concatenated_merged_abundance_blast_taxa.txt
#             fi

#             if [ ! -s "/home/ubuntu/myMetagenomics/Trial/04_abundance/class_family_genus_species.txt" ]; then
#                 touch /home/ubuntu/myMetagenomics/Trial/04_abundance/class_family_genus_species.txt
#             fi
#         fi
#         """

# #所有samples的abundance计算完后再执行
# rule merge_virus_abundance_taxa_host:
#     input:
#         merged_abundance_blast_taxa = join(abandunce_dir, "all_concatenated_merged_abundance_blast_taxa.txt"),
#         class_family_genus_species_file = join(abandunce_dir, "class_family_genus_species.txt"),
#     output:
#         virus_host_dictionary_with_common = join(abandunce_dir, "my_virus_host_dictionary_with_common.csv"),
#         all_concatenated_file = join(abandunce_dir, "quantification_result.txt"),
#     conda: "base"
#     params:
#         file_dir = abandunce_dir,
#         virushostdb_file = '/home/ubuntu/ReferenceData/virushostdb.tsv',

#         #这两个需要自己时常更新，因为有些病毒信息是随着样本新出现的，需要自己补充(如果发现有NA的需要补充完后再跑一次)
#         my_virus_host_dictionary = '/home/ubuntu/ReferenceData/my_virus_host_dictionary.csv',
#         host_dict_file = '/home/ubuntu/ReferenceData/host_dictionarty.tsv',
#     shell:
#         """ 
#         #combine all merged_abundance_blast_taxa_grouped_by_species.txt
#         bash /home/ubuntu/Tools/Scripts/combine_all_merged.sh {params.file_dir} {input.merged_abundance_blast_taxa}

#         #generate species_host file
#         python /home/ubuntu/Tools/Scripts/generate_species_host.py {input.class_family_genus_species_file} {params.virushostdb_file} {params.my_virus_host_dictionary}

#         ##To add a new "host_common" column to your my_virus_host_dictionary.csv
#         python /home/ubuntu/Tools/Scripts/host_categorize.py {params.my_virus_host_dictionary} {params.host_dict_file} {output.virus_host_dictionary_with_common}

#         #merge host_common and all_concatenated_merged_abundance_blast_taxa.txt
#         python /home/ubuntu/Tools/Scripts/combine_virus_abundance_host.py {input.merged_abundance_blast_taxa} {output.virus_host_dictionary_with_common} {output.all_concatenated_file}

#         grep -viE "phage|bacteriophage|Caudovir|Microviridae|Microvirus|Myoviridae|Prokaryotic|Siphoviridae|Siphovirus|Podoviridae|Podovirus|Inoviridae|leviviridae|Herelleviridae|Ackermannviridae|Crassvirales" {output.all_concatenated_file} > {output.all_concatenated_file}.no.phages

#         """

# # 这里为病毒添加DNA还是RNA的注释信息
# rule merge_dna_rna_annotation:
#     input:
#         all_concatenated_file_no_phage = join(abandunce_dir, "quantification_result.txt.no.phages"),
#     output:
#         all_concatenated_file_no_phage_with_genome_annotation = join(abandunce_dir, "quantification_result.txt.no.phages_with_genome.txt"),
#     conda: "base"
#     params:
#         ictv_master_species = "/home/ubuntu/ReferenceData/ICTV_MS_39_v4.tsv",
#     shell:
#         """
#         python /home/ubuntu/Tools/Scripts/merge_dna_rna_annotation.py {input.all_concatenated_file_no_phage} {params.ictv_master_species} {output.all_concatenated_file_no_phage_with_genome_annotation}
#         """

# #在上一步中，我已经把相同Species的病毒合并了（RPM, Mapped_reads），所以这一步不需要再合并
# #Filter on read count of each virus18. If the total read count of a specific virus in a specific library of is less than 0.1% the highest read count for that virus within the same sequencing lane, then it is considered as a false positive due to index-hopping. 
# #这里过滤掉了几个被RT-PCR检测出了流感的病毒，因此暂时不用这些Alphainfluenzavirus influenzae
# #仍然用quantification_result.txt.no.phages_with_genome.txt作为输入
# rule remove_index_hopping:
#     input:
#         all_concatenated_file = join(abandunce_dir, "quantification_result.txt.no.phages_with_genome.txt"),
#         samples_metadata = join("/home/ubuntu/myMetagenomics/Trial", "Metagenomic_Samples_Annotation_file_with_groups.csv")
#     output:
#         all_concatenated_file_no_index_hopping = join(abandunce_dir, "quantification_result.txt.no.phages_with_genome.txt.indexhopping"),
#         filtered_index_hopping_file = join(abandunce_dir, "filtered.index.hopping.txt")
#     conda: "base"
#     shell:
#         """
#         python /home/ubuntu/Tools/Scripts/filter_index_hopping.py {input.all_concatenated_file} {input.samples_metadata} {output.all_concatenated_file_no_index_hopping} {output.filtered_index_hopping_file}
#         """