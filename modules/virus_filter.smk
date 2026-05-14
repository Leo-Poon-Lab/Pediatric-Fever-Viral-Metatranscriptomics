rule reads_contigs_integration:
    input:
        reads_virus_summary = join(reads_taxa_kraken_kaiju_dir, "virus_summary.tsv"),
        contigs_virus_summary = join(virus_identification_dir,"all_combined_virus_results_lineages.tsv")
    output:
        reads_contgis_virus_summary = join(reads_contigs_dir,"virus_summary_with_contigs.tsv"),
    conda:
        "envs/workflow.yaml"
    shell:
        """
            #4. Reads Contigs Info
            python scripts/python/combine_reads_contigs_info.py {abandunce_dir} {input.reads_virus_summary} {input.contigs_virus_summary} {output.reads_contgis_virus_summary}

        """

rule filter:
    input:
        reads_contgis_virus_summary = join(reads_contigs_dir,"virus_summary_with_contigs.tsv"),
    output:
        reads_contgis_virus_summary_index_hopping = join(reads_contigs_dir,"virus_summary_with_contigs_indexhopping.tsv"),
        reads_contgis_virus_summary_index_hopping_filtered_1 = join(reads_contigs_dir,"virus_summary_with_contigs_indexhopping.filtered.1.tsv"),
        reads_contgis_virus_summary_index_hopping_filtered_final = join(reads_contigs_dir,"virus_summary_with_contigs_indexhopping.filtered.final.tsv"),
        reads_contgis_virus_summary_index_hopping_filtered_final_dnarna = join(reads_contigs_dir,"virus_summary_with_contigs_indexhopping.filtered.final.dnarna.tsv")
    params:
        patient_profile = "config/patient_metadata.tsv",
        reads_contgis_virus_summary_index_hopping_removed = join(reads_contigs_dir, "virus_summary_with_contigs_indexhopping_removed.tsv"),
        contamination_annotation_file = Path(config["databases"]["annotation"]["contamination"]),
        ICTV_MS_39_v4 = Path(config["databases"]["taxonomy"]["ictv"]),
    conda:
        "envs/workflow.yaml"
    shell:
        """
            # Final reads support and RPM
            # Index-Hopping 0.1% threhold for  0.1% of the highest read count for that virus within the same sequencing lane
            python scripts/python/filter_index_hopping_reads.py {input.reads_contgis_virus_summary} {params.patient_profile} {output.reads_contgis_virus_summary_index_hopping} {params.reads_contgis_virus_summary_index_hopping_removed}

            # Reagent contamination; Phages; RPM >= 1; Number_of_supporting_reads >= 3; Average_identity >= 90; Average alignment length >= 70;
            # 对于有Contigs支持的taxon，要求至少有 contig_length >= 300 and mapped_reads >= 3 and rpm_mapping_contigs >= 1
            python scripts/python/comprehensive_filter.py {output.reads_contgis_virus_summary_index_hopping} {output.reads_contgis_virus_summary_index_hopping_filtered_1}

            #Contamination reference
            python scripts/python/filter_contaminants.py {output.reads_contgis_virus_summary_index_hopping_filtered_1} {params.contamination_annotation_file} {output.reads_contgis_virus_summary_index_hopping_filtered_final}

        
            #merge_dna_rna_annotation
            python  scripts/python/merge_dna_rna_annotation.py {output.reads_contgis_virus_summary_index_hopping_filtered_final} {params.ICTV_MS_39_v4} {output.reads_contgis_virus_summary_index_hopping_filtered_final_dnarna}

        """