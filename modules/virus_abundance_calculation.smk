checkpoint contigs_alignment:
    input:
        rmhost_r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
        rmhost_r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz"),
        virus_contigs = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs.fa")
    output:
        aligned_reads = join(abandunce_dir, "{sample}/virus_aligned.bam"),
    conda: "../envs/workflow.yaml"
    threads: 5
    params:
        index_dir = join(abandunce_dir, "{sample}/viral_contigs_index"),
        index = join(abandunce_dir, "{sample}/viral_contigs_index/viral_contigs_index"),
        sensitivity = "--very-sensitive"
    shell:
        """
        if [ ! -s {input.virus_contigs} ]; then
            touch {output.aligned_reads}
            exit 0
        else
            mkdir -p {params.index_dir}
            bowtie2-build {input.virus_contigs} {params.index}
            bowtie2 -p {threads} -x {params.index} -1 {input.rmhost_r1} -2 {input.rmhost_r2} \
            | samtools sort -O bam -@ {threads} -m 10G -o - > {output.aligned_reads}
        fi
       
        """

rule coverage_stat:
    input:
        aligned_reads = join(abandunce_dir, "{sample}/virus_aligned.bam")
    output:
        coverage_log = join(abandunce_dir, "{sample}/weesam_virus_coverage.txt"),
        coverage_sorted = join(abandunce_dir, "{sample}/coverage_sorted.tsv"),
    conda: "../envs/workflow.yaml"
    threads: 10
    params:
        file_dir = join(abandunce_dir, "{sample}"),
        wee_sam_html = join(abandunce_dir, "{sample}/wee_sam.html"),
    shell:
        """
        if [ ! -s {input.aligned_reads} ]; then
            touch {output.coverage_log}
            touch {output.coverage_sorted}
            exit 0
        else        
            /home/ubuntu/Tools/weeSAM/weeSAM-master/weeSAM --bam {input.aligned_reads} \
            --out {output.coverage_log} --html {params.wee_sam_html}  --overwrite  --pretty > {params.file_dir}/weesam.log 2>&1

            cat {output.coverage_log} | cut -f1,2,3,4,5,8 | awk -v FS="\t" -v OFS="\t" '{{if($5>=0) print $0}}' | sort -t $'\t' -k5nr -k6nr > {params.file_dir}/coverage_sorted.tsv
        fi

        """

rule virus_abundance_estimation:
    input:
        aligned_reads = join(abandunce_dir, "{sample}/virus_aligned.bam"),
        virus_coverage_file = join(abandunce_dir, "{sample}/coverage_sorted.tsv"),
    output:
        virus_abundance = join(abandunce_dir, "{sample}/virus_abundance.txt"),
    conda: "../envs/workflow.yaml"
    params:
        file_dir = join(abandunce_dir, "{sample}"),
        sample_name = "{sample}"
    shell:
        """
        total_reads=$(samtools view -c {input.aligned_reads})
        echo -e "Contig\tLength\tMapped_Reads\tBreadth\tCoverage\tAvg_Depth\tRPM" > {output.virus_abundance}

        awk -v total_reads=$total_reads \
        'BEGIN {{ OFS="\t" }} {{ rpm = ($3 / total_reads) * 1000000; \
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%.2f\\n", $1, $2, $3, $4, $5, $6, rpm) }}' \
        {input.virus_coverage_file} >> {output.virus_abundance}
        
        """

rule merge_virus_abundance_taxa:
    input:
        virus_abundance = join(abandunce_dir, "{sample}/virus_abundance.txt"),
        virus_matches = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs_all_no_phage.tsv"),
    output:
        virus_abundance_taxa = join(abandunce_dir, "{sample}/merged_abundance_blast_taxa.txt"),
        virus_abundance_taxa_grouped = join(abandunce_dir, "{sample}/merged_abundance_blast_taxa_grouped_by_species.txt"),
    conda: "../envs/workflow.yaml"
    params:
        sample_name = "{sample}"
    shell:
        """ 
        if [ ! -s {input.virus_matches} ]; then
            touch {output.virus_abundance_taxa}
            touch {output.virus_abundance_taxa_grouped}
            exit 0
        else  
            python scripts/python/contigs/merge_abundance_blast_taxa.py {params.sample_name} {input.virus_abundance} {input.virus_matches} {output.virus_abundance_taxa}
            
            python scripts/python/contigs/sum_rpm_and_group_species.py {output.virus_abundance_taxa} {output.virus_abundance_taxa_grouped}

        fi
        """

#所有samples的abundance计算完后再执行
rule merge_virus_abundance_taxa_all:
    input:
        merged_abundance_blast_taxa = expand(
            join(abandunce_dir, "{sample}/merged_abundance_blast_taxa_grouped_by_species.txt"),
            sample=SAMPLES
        )
    output:
        merged_abundance_blast_taxa = join(abandunce_dir, "all_concatenated_merged_abundance_blast_taxa.txt"),
        merged_abundance_blast_taxa_no_phage = join(abandunce_dir, "all_concatenated_merged_abundance_blast_taxa.no.phages.txt")
    conda: "../envs/workflow.yaml"
    shell:
        """ 
        #combine all merged_abundance_blast_taxa_grouped_by_species.txt
        bash scripts/shell/contigs/combine_all_merged.sh {abandunce_dir} {output.merged_abundance_blast_taxa}

        if [ ! -s {output.merged_abundance_blast_taxa} ]; then
            touch {output.merged_abundance_blast_taxa}
            touch {output.merged_abundance_blast_taxa_no_phage}
        else
            grep -viE "phage|bacteriophage|Caudovir|Microviridae|Microvirus|Myoviridae|Prokaryotic|Siphoviridae|Siphovirus|Podoviridae|Podovirus|Inoviridae|leviviridae|Herelleviridae|Ackermannviridae|Crassvirales" {output.merged_abundance_blast_taxa} > {output.merged_abundance_blast_taxa_no_phage}

        fi

        """

rule merge_dna_rna_annotation:
    input:
        merged_abundance_blast_taxa_no_phage = join(abandunce_dir, "all_concatenated_merged_abundance_blast_taxa.no.phages.txt")
    output:
        merged_abundance_blast_taxa_no_phage_with_genome_annotation = join(abandunce_dir, "all_concatenated_merged_abundance_blast_taxa.no.phages_with_genome.txt"),
    conda: "../envs/workflow.yaml"
    params:
        ictv_master_species = Path(config["databases"]["taxonomy"]["ictv"]),
    shell:
        """
        python scripts/python/contigs/merge_dna_rna_annotation.py {input.merged_abundance_blast_taxa_no_phage} {params.ictv_master_species} {output.merged_abundance_blast_taxa_no_phage_with_genome_annotation}
        """

#Filter on read count of each virus18. If the total read count of a specific virus in a specific library of is less than 0.1% the highest read count for that virus within the same sequencing lane, then it is considered as a false positive due to index-hopping. 
rule remove_index_hopping:
    input:
        merged_abundance_blast_taxa_no_phage_with_genome_annotation = join(abandunce_dir, "all_concatenated_merged_abundance_blast_taxa.no.phages_with_genome.txt"),
        samples_metadata = join("config/", "patient_metadata.tsv")
    output:
        merged_abundance_blast_taxa_no_phage_with_genome_annotation_no_index_hopping = join(abandunce_dir, "all_concatenated_merged_abundance_blast_taxa.no.phages_with_genome.txt.indexhopping"),
        filtered_index_hopping_file = join(abandunce_dir, "filtered.index.hopping.txt")
    conda: "base"
    shell:
        """
        python scripts/python/contigs/filter_index_hopping.py {input.merged_abundance_blast_taxa_no_phage_with_genome_annotation} {input.samples_metadata} {output.merged_abundance_blast_taxa_no_phage_with_genome_annotation_no_index_hopping} {output.filtered_index_hopping_file}
        """