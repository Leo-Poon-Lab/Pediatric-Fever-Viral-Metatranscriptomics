#所有samples的abundance计算完后再执行,这一步进行了筛选和过滤
rule filter_RPM_contig_length:
    input:
        all_concatenated_file_no_index_hopping = join(abandunce_dir, "quantification_result.txt.no.phages"),
    output:
        all_concatenated_file_no_index_hopping_rpm_length = join(virus_component_analysis_dir, "quantification_result.txt.no.phages.RPM_length"),
        all_concatenated_file_no_index_hopping_rpm_length_removed = join(virus_component_analysis_dir, "quantification_result.txt.no.phages.RPM_length_removed"),
    conda: "base"
    params:
        RPM_threshold = 1,
        contig_length = 300,
        file_dir = virus_component_analysis_dir, 
    shell:
        """ 
        # 这个脚本可以调节RPM的筛选阈值
        bash /home/ubuntu/Tools/Scripts/filter_RPM_contig_length.sh {input.all_concatenated_file_no_index_hopping} {params.RPM_threshold} {params.contig_length} {output.all_concatenated_file_no_index_hopping_rpm_length} {output.all_concatenated_file_no_index_hopping_rpm_length_removed}        

        """
#这一步是将RT_PCR的结果加入到我们的结果中，作为筛选contigs的金标准
rule add_RT_PCR_to_blast_RPM:
    input:
        all_combined_diamond_blastn_without_phage = join(virus_identification_dir, "all_combined_virus_results.no.phages.txt"),
        all_concatenated_file_no_phage_with_genome_annotation = join(abandunce_dir, "quantification_result.txt.no.phages_with_genome.txt")
    output:
        all_combined_diamond_blastn_without_phage_with_RT_PCR = join(virus_component_analysis_dir, "all_combined_virus_results.no.phages.with_RT_PCR.txt"),
        all_concatenated_file_no_phage_with_genome_annotation_with_RT_PCR = join(virus_component_analysis_dir, "quantification_result.txt.no.phages_with_genome.with_RT_PCR.txt"),
    params:
        patient_metadata = "/home/ubuntu/myMetagenomics/Trial/data/Patient_Profile_Specimen_new_sequencing.txt",
    conda: "base"
    shell:
        """ 
        python /home/ubuntu/Tools/Scripts/merge_virus_results_RT-PCR.py {input.all_combined_diamond_blastn_without_phage} {params.patient_metadata} {output.all_combined_diamond_blastn_without_phage_with_RT_PCR}
        python /home/ubuntu/Tools/Scripts/merge_virus_results_RT-PCR.py {input.all_concatenated_file_no_phage_with_genome_annotation} {params.patient_metadata} {output.all_concatenated_file_no_phage_with_genome_annotation_with_RT_PCR}

        """
#将Kraken2的reads分析结果以及
rule combine_kraken_reads_contigs_blast_results:
    input:
        combined_sample_virus_confirmed_noPhages = join(reads_taxa_kranken2_dir, "all_virus_confirmed_noPhages.tsv"),
        all_concatenated_file_no_phage_with_genome_annotation_with_RT_PCR = join(virus_component_analysis_dir, "quantification_result.txt.no.phages_with_genome.with_RT_PCR.txt"),
    output:
        all_combined_kraken_diamond_results = join(virus_component_analysis_dir, "all_combined_kraken_diamond_results.txt"),
    params:
        patient_metadata = "/home/ubuntu/myMetagenomics/Trial/data/Patient_Profile_Specimen_new_sequencing.txt",
    conda: "base"
    shell:
        """
        python /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/combine_kraken_reads_blast_contigs_tables.py

        """

#将rpm和blast的所有结果联合分析，目前需要手动挑选blast结果
rule combine_blast_rpm:
    input:
        all_concatenated_file_no_index_hopping_rpm_length = join(virus_component_analysis_dir, "quantification_result.txt.no.phages.RPM_length"),
        all_combined_diamond_blastn_without_phage = join(virus_identification_dir, "all_combined_virus_results.no.phages.txt"),
    output:
        all_combined_file = join(virus_component_analysis_dir, "quantification_rpm_blast_results.tsv"),
    conda: "base"
    shell:
        """ 
        python /home/ubuntu/Tools/Scripts/combine_blast_rpm.py {input.all_concatenated_file_no_index_hopping_rpm_length} {input.all_combined_diamond_blastn_without_phage} {output.all_combined_file}

        """

#这一步对可能存在的false positive (主要是那篇文献提到的可能的污染物列表)进行最终地过滤
rule filter_false_positives:
    input:
        all_combined_file = join(virus_component_analysis_dir, "quantification_rpm_blast_results.tsv" ),
        all_concatenated_file_no_index_hopping_rpm_length = join(virus_component_analysis_dir, "quantification_result.txt.no.phages.RPM_length"),
    output:
        output_contam_file = join(virus_component_analysis_dir, "potential_contamination_contigs.txt"),
        output_results_file = join(virus_component_analysis_dir, "quantification_rpm_blast_rusults.tsv_con_filtered"),
        all_concatenated_file_no_index_hopping_rpm_length_filtered = join(virus_component_analysis_dir, "quantification_result.txt.no.phages.RPM_length_filtered"),
    conda: "base"
    params:
        file_dir = virus_component_analysis_dir,
        contam_file = "/home/ubuntu/myMetagenomics/Trial/contamination_lab_known_taxID.txt",
    shell:
        """ 
        #对contigs按照一定条件进行过滤
        python /home/ubuntu/Tools/Scripts/filter_potential_false_posititives_contigs.py {params.contam_file} {input.all_combined_file} {output.output_contam_file} {output.output_results_file}

        #根据已知的contamination进行过滤
        python /home/ubuntu/Tools/Scripts/filter_known_contamination.py {params.contam_file} {input.all_combined_file} {output.output_contam_file} {output.output_results_file}

        #将contamination筛选后的Contig从文件quantification中提取出来用于后面的heatmap
        python /home/ubuntu/Tools/Scripts/filter_quantification_results.py {output.output_results_file} {input.all_concatenated_file_no_index_hopping_rpm_length} {output.all_concatenated_file_no_index_hopping_rpm_length_filtered}

        sed 's/,/\t/g' {output.all_concatenated_file_no_index_hopping_rpm_length_filtered} > {output.all_concatenated_file_no_index_hopping_rpm_length_filtered}.tsv
        sed 's/,/\t/g' {output.output_results_file} > {output.output_results_file}.tsv
        sed 's/,/\t/g' {output.output_contam_file} > {output.output_contam_file}.tsv

        """

#上一步的结果还需再看一下，最好手动再检查一下，最后再进行heatmap
#所有samples的abundance计算完后再执行
rule heatmap:
    input:
        all_combined_file_filtered = join(virus_component_analysis_dir, "quantification_result.txt.no.phages.RPM_length_filtered"),
    output:
        all_concatenated_file = join(virus_component_analysis_dir, "OK.txt"),
    conda: "R"
    params:
        file_dir = virus_component_analysis_dir,
        all_domain_host_meta = join(virus_component_analysis_dir, "all_domain_host_meta.csv"),
        all_table_abundance = join(virus_component_analysis_dir, "all_table_abundance.txt"),
        top_metadata = join(virus_component_analysis_dir, "top_metadata.csv"),
    shell:
        """ 
        bash  /home/ubuntu/Tools/Scripts/graphics/heat/snakemake/generate_needed_files.sh {input.all_combined_file_filtered} {params.all_domain_host_meta} {params.all_table_abundance} {params.top_metadata}
        Rscript /home/ubuntu/Tools/Scripts/graphics/heat/snakemake/draw_heatmap_virus_no_order_batch.R

        """
# ★★★
#在reads分析完后，直接来这里
rule heatmap_new_trial:
    input:
        all_combined_file_filtered = join("/home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/new_trial/", "virus_summary_with_contigs_indexhopping.3.filtered.added.dna-rna.contigs.tsv"),
    output:
        all_concatenated_file = join(virus_component_analysis_dir, "OK.txt"),
    conda: "R"
    params:
        file_dir = virus_component_analysis_dir,
        all_domain_host_meta = join(virus_component_analysis_dir, "all_domain_host_meta.csv"),
        all_table_abundance = join(virus_component_analysis_dir, "all_table_abundance.txt"),
        top_metadata = join(virus_component_analysis_dir, "top_metadata.csv"),
    shell:
        """ 
        # 这里做所有病毒的heatmap
        cd /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/new_trial
        bash  /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/new_trial/generate_needed_files_new_trial.sh
        Rscript /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/new_trial/draw_heatmap_virus_new_trial.R

        Rscript draw_heatmap_virus_new_trial_genus.R 
        # 这里做所有病毒的heatmap，横向展示
        Rscript draw_heatmap_virus_new_trial_horizontal.R 

        ###### Focus RT-PCR的样本结果, 画上面热图的子集 ########
        cd /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/new_trial/RT-PCR
        bash  /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/new_trial/RT-PCR/generate_needed_files_new_trial_RT-PCR.sh
        Rscript /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/new_trial/RT-PCR/draw_heatmap_virus_new_trial_RT_PCR.R

        # 这里专注于RT-PCR那一系列呼吸道病毒的检测对比，RT-PCR的标记*
        # Rscript /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/new_trial/RT-PCR/draw_RT_PCR_vs_metagenomics.R
        python /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/new_trial/RT-PCR/collect_RTPCR_respirotory_viruses.py 
        Rscript /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/new_trial/RT-PCR/draw_RT_PCR_vs_metagenomics.R

        """
# 新的分析
rule heatmap_new_trial_virus:
    shell:
        """ 

   
        """


#检测RT-PCR与Reads-level和Contigs-level的检出率
rule rt_pcr_stat:
    input:
        patient_file = "/home/ubuntu/myMetagenomics/Trial/data/Patient_Profile_Specimen_new.txt",
        reads_results = join(reads_taxa_kraken_kaiju_dir, "{sample}/all_kraken_kaiju_combined_report_blastn_PCR.tsv"),
        contigs_results = join(reads_taxa_kraken_kaiju_dir, "{sample}/all_kraken_kaiju_combined_report_blastn_PCR.tsv"),
    output:
        stat_results = join(virus_component_analysis_dir, "RT-PCR/all_table_abundance.txt"),
    conda: "base"
    params:
        file_dir = virus_component_analysis_dir,
        contam_file = "/home/ubuntu/myMetagenomics/Trial/contamination_lab_known_taxID.txt",
    shell:
        """ 
        #对contigs按照一定条件进行过滤
        python /home/ubuntu/myMetagenomics/Trial/05_virus_component_analysis/RT-PCR/calc_sensitive_methods.py {input.patient_file} {input.reads_results} {input.contigs_results}

        """
