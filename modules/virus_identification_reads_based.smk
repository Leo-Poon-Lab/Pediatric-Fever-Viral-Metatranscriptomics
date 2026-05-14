rule kraken2:
    input:
        r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
        r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz")
       
    output:
        kreport = join(reads_taxa_kranken2_dir,"{sample}/{sample}.kreports"),
        output_file = join(reads_taxa_kranken2_dir,"{sample}/{sample}.output.txt"),
    threads: 10
    params:
        database = Path(config["databases"]["taxonomy"]["kraken2uniq"]["db"]),
        minimum_hit = 3,
        confidence = 0.05
    conda:
        "envs/workflow.yaml"
    shell:
        """
        # --report-minimizer-data加这个参数就是KrakenUniq
        kraken2 --db {params.database} --paired --threads {threads} --report {output.kreport} --output {output.output_file} \
        --report-minimizer-data --minimum-hit-groups {params.minimum_hit} \
        --confidence {params.confidence} \
        {input.r1} {input.r2}
        """

rule bracken:
    input:
        kreport = join(reads_taxa_kranken2_dir,"{sample}/{sample}.kreports"),
    output:
        breport =  join(reads_taxa_kranken2_dir,"{sample}/{sample}.breports"),
        output_bracken = join(reads_taxa_kranken2_dir,"{sample}/{sample}.bracken"),
    params:
        database = Path(config["databases"]["taxonomy"]["kraken2uniq"]["db"]),
        read_len = 100,
        abundance_level = "S"
    conda:
        "envs/workflow.yaml"
    shell:
        "bracken -d {params.database} -i {input.kreport} -w {output.breport} -o {output.output_bracken} -r {params.read_len} -l {params.abundance_level}"

rule krona_plots:
    input:
        breport = join(reads_taxa_kranken2_dir,"{sample}/{sample}.breports"),
    output:
        output_krona = join(reads_taxa_kranken2_dir,"{sample}/{sample}.b.krona.txt"),
    log:
        html =  join(reads_taxa_kranken2_dir,"{sample}/{sample}.krona.html"),
    conda:
        "envs/workflow.yaml"
    shell:
        """
        python scripts/tools/KrakenTools-1.2/kreport2krona.py --no-intermediate-ranks -r {input.breport} -o {output.output_krona}
        ktImportText {output.output_krona} -o {log}
        """

rule kaiju:
    input:
        r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
        r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz")   
    output:
        kaiju_out = join(reads_taxa_kaiju_dir, "{sample}/kaiju_greedy.out"),
        kaiju_report = join(reads_taxa_kaiju_dir, "{sample}/kaiju_greedy_report.tsv")
    threads: 10
    params:
        nodes = Path(config["databases"]["taxonomy"]["kaiju"]["nodes"]),
        names = Path(config["databases"]["taxonomy"]["kaiju"]["names"]),
        db = Path(config["databases"]["taxonomy"]["kaiju"]["db"]),  # specify your Kaiju database path

        nodes_refseq = Path(config["databases"]["taxonomy"]["kaiju"]["refseq_nodes"]),
        names_refseq = Path(config["databases"]["taxonomy"]["kaiju"]["refseq_names"]),
        db_refseq = Path(config["databases"]["taxonomy"]["kaiju"]["refseq_db"]),  # specify your Kaiju database path
        mismatch = 5,  # Number of mismatches allowed in Greedy mode (default: 3)
        evalue = 0.01  # Minimum E-value in Greedy mode
    conda:
        "envs/workflow.yaml"
    shell:
        """
        # Run Kaiju classification with sensitive parameters
        # Greedy run mode yields a higher sensitivity compared with MEM mode
        kaiju -t {params.nodes} -f {params.db} -i {input.r1} -j {input.r2} -o {output.kaiju_out} -z {threads} \
              -a greedy -e {params.mismatch}

        # Generate taxonomy report
        kaiju2table -t {params.nodes} \
                    -n {params.names} \
                    -r species \
                    -l superkingdom,class,family,genus,species -e \
                    -o {output.kaiju_report} \
                    {output.kaiju_out}

        """

rule combine_kaiju_kraken:
    input:
        kraken_report = join(reads_taxa_kranken2_dir,"{sample}/{sample}.kreports"),
        kaiju_report = join(reads_taxa_kaiju_dir, "{sample}/kaiju_greedy_report.tsv")
    output:
        kraken_report_reformatted_virus = join(reads_taxa_kranken2_dir,"{sample}/{sample}.kreports.reformatted.virus"),
        Kraken_Kaiju_combined_report = join(reads_taxa_kraken_kaiju_dir,"{sample}/kraken_kaiju_combined_report.tsv")
    params:
        work_dir = join(reads_taxa_kraken_kaiju_dir, "{sample}"),
    conda:
        "envs/workflow.yaml"
    shell:
        """
        #Kraken
        # 将Kraken的层级输出改为Kaiju的层级输出，易于理解以及后续合并分析
        awk '$6 != "F1" && $6 != "F2" && $6 != "F3" && $6 != "G1" && $6 != "G2" && $6 != "G3" && $6 != "S1" && $6 != "S2" && $6 != "S3" {{ print }}' {input.kraken_report} > {input.kraken_report}.new
        python scripts/tools/KrakenTools-1.2/reformat_kraken_kreports_virus.py {input.kraken_report}.new {output.kraken_report_reformatted_virus}
        
        # Ensure reformatted virus file exists
        if [ ! -f {output.kraken_report_reformatted_virus} ]; then
            touch {output.kraken_report_reformatted_virus}
        fi
        
        # Filter phages from Kraken results
        grep -i "virus" {output.kraken_report_reformatted_virus} | grep -vi -e "phage" -e "Caudoviricete" -e "Microviridae" -e "Microvirus" -e "Myoviridae" -e "Prokaryotic" -e "Siphoviridae" -e "Siphovirus" -e "Podoviridae" -e "Podovirus" -e "Inoviridae" -e "leviviridae" -e "Herelleviridae" -e "Ackermannviridae" \
            > {output.kraken_report_reformatted_virus}.noPhages

        # Kaiju processing
        grep -i "virus" {input.kaiju_report} | grep -vi -e "phage" -e "Caudoviricete" -e "Microviridae" -e "Microvirus" -e "Myoviridae" -e "Prokaryotic" -e "Siphoviridae" -e "Siphovirus" -e "Podoviridae" -e "Podovirus" -e "Inoviridae" -e "leviviridae" -e "Herelleviridae" -e "Ackermannviridae" \
            > {input.kaiju_report}.virus.noPhages
        
        mkdir -p {params.work_dir}
        # Handle empty results
        if [ ! -s {output.kraken_report_reformatted_virus}.noPhages ] && [ ! -s {input.kaiju_report}.virus.noPhages ]; then
            echo "No viral entries found. Creating empty output."
            touch {output.Kraken_Kaiju_combined_report}
        else
            # combine kaiju and kranken results
            # Union Set
            python scripts/python/reads/combine_kaiju_kraken.py \
                {output.kraken_report_reformatted_virus}.noPhages \
                {input.kaiju_report}.virus.noPhages \
                {input.kraken_report} \
                {output.Kraken_Kaiju_combined_report}
        fi
        
        """

rule extract_kaiju_kraken_reads:
    input:
        combined_report = join(reads_taxa_kraken_kaiju_dir,"{sample}/kraken_kaiju_combined_report.tsv"),
        r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
        r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz")
    output:
        reads_extracted_dir = directory(join(reads_taxa_kraken_kaiju_dir, "{sample}/extracted_reads"))
    params:
        kraken_output = join(reads_taxa_kranken2_dir,"{sample}/{sample}.output.txt"),
        kreport = join(reads_taxa_kranken2_dir,"{sample}/{sample}.kreports"),
        kaiju_out = join(reads_taxa_kaiju_dir, "{sample}/kaiju_greedy.out"),
    conda:
        "envs/workflow.yaml"
    shell:
        """
        python scripts/python/reads/extract_reads_from_kraken_kaiju_combined_report.py \
        --combined-report {input.combined_report} \
        --input-r1 {input.r1} --input-r2 {input.r2} \
        --kraken-output {params.kraken_output} \
        --kraken-report {params.kreport} \
        --kaiju-output {params.kaiju_out}  \
        --output-dir {output.reads_extracted_dir} --deduplicate

        """

rule combined_reads_BLASTn_confirmed:
    input:
        reads_extracted_dir = join(reads_taxa_kraken_kaiju_dir, "{sample}/extracted_reads")
    output:
        blast_results_dir = directory(join(reads_taxa_kraken_kaiju_dir, "{sample}/blast_results"))
    params:
        reads_renamed_dir = join(reads_taxa_kraken_kaiju_dir, "{sample}/renamed_extracted_reads"),

        database = Path(config["databases"]["nucleotide"]["nt"]),
        task = "blastn",
        evalue = 1e-5,
        threads = 5,
    conda:
        "envs/workflow.yaml"
    shell:
        """
        # Validate input directory exists
        if [ ! -d "{input.reads_extracted_dir}" ]; then
            echo "Error: Input directory {input.reads_extracted_dir} does not exist" >&2
            exit 1
        fi

        # Create output directories if they don't exist
        mkdir -p "{params.reads_renamed_dir}" || true
        mkdir -p "{output.blast_results_dir}" || true

        # Rename reads header with existence check
        if [ -n "$(ls -A {input.reads_extracted_dir}/*.fa 2>/dev/null)" ]; then
            python scripts/python/reads/rename_and_combine_kraken_kaiju_reads.py \
                --input_dir "{input.reads_extracted_dir}" \
                --output_dir "{params.reads_renamed_dir}"
        else
            echo "Warning: No FASTA files found in {input.reads_extracted_dir}"
            touch "{params.reads_renamed_dir}/.empty"
        fi

        # Run BLASTn validation if renamed files exist
        if [ -n "$(ls -A {params.reads_renamed_dir}/*.fa 2>/dev/null)" ]; then
            python scripts/python/reads/blastn_runner.py \
                --input_dir "{params.reads_renamed_dir}" \
                --output_dir "{output.blast_results_dir}" \
                --database "{params.database}" \
                --task "{params.task}" \
                --evalue "{params.evalue}" \
                --threads "{params.threads}"  
        else
            echo "Warning: No renamed FASTA files found for BLASTn"
            touch "{output.blast_results_dir}/.empty"
        fi
        """

rule integration_analsis:
    input:
        blast_results_dir = join(reads_taxa_kraken_kaiju_dir, "{sample}/blast_results"),
        Kraken_Kaiju_combined_report = join(reads_taxa_kraken_kaiju_dir,"{sample}/kraken_kaiju_combined_report.tsv")
    output:
        Kraken_Kaiju_combined_report_class = join(reads_taxa_kraken_kaiju_dir,"{sample}/kraken_kaiju_combined_report_with_class.tsv"),
        Kraken_Kaiju_combined_report_class_blast = join(reads_taxa_kraken_kaiju_dir,"{sample}/kraken_kaiju_combined_report_with_class_blastn.tsv"),
    params:
        temp_file= join(reads_taxa_kraken_kaiju_dir, "{sample}/temp_taxids.txt"),
        temp_taxonkit_output= join(reads_taxa_kraken_kaiju_dir, "{sample}/temp_taxonkit_output.txt"),
        temp_all_unique_taxids = join(reads_taxa_kraken_kaiju_dir, "{sample}/all_unique_taxids.txt")
    conda:
        "base"
    shell:
        """
        # Validate input directory exists
        if [ ! -d "{input.blast_results_dir}" ]; then
            echo "Error: Input directory {input.blast_results_dir} does not exist" >&2
            exit 1
        fi

        # ★★★
        # 主体思路是对照query和subject的Viral Class, TaxID部分看是否一致
        # Combine BLAST results, add lineage info for Subject and Query
        bash scripts/shell/reads/add_lineage_to_blast.sh {input.blast_results_dir}
        # 添加Class信息
        bash scripts/shell/reads/add_lineage_class_to_Kraken_Kaiju_combined_report.sh {input.Kraken_Kaiju_combined_report} {output.Kraken_Kaiju_combined_report_class} {params.temp_file} {params.temp_taxonkit_output} {params.temp_all_unique_taxids}

        # Add if BLASTn supports Kraken or Kaiju classification
        python scripts/python/reads/add_blast_confirmed_column.py --blast_dir {input.blast_results_dir} \
        --combined_report {output.Kraken_Kaiju_combined_report_class} \
        --output_report {output.Kraken_Kaiju_combined_report_class_blast} \
        --evalue_threshold 1e-5 \
        --identity_threshold 70 

        """

rule integration_analsis_2:
    input:
        Kraken_Kaiju_combined_report_class_blast_results = expand(
            join(reads_taxa_kraken_kaiju_dir, "{sample}/kraken_kaiju_combined_report_with_class_blastn.tsv"),
            sample=SAMPLES
        )
    output:
        all_Kraken_Kaiju_combined_blast_results = join(reads_taxa_kraken_kaiju_dir, "all_kraken_kaiju_combined_report_class_blastn.tsv"),
        all_blast_results_virus_lineages = join(reads_taxa_kraken_kaiju_dir, "all_blast_results_virus_lineages.tsv"),
        virus_summary = join(reads_taxa_kraken_kaiju_dir, "virus_summary.tsv")
    params:
        all_Kraken_Kaiju_combined_blast_results_pcr = join(reads_taxa_kraken_kaiju_dir, "all_kraken_kaiju_combined_report_class_blastn_PCR.tsv"),
    conda:
        "base"
    shell:
        """
            find {reads_taxa_kraken_kaiju_dir} -name kraken_kaiju_combined_report_with_class_blastn.tsv | xargs awk "FNR==1 && NR!=1 {{next;}}{{print}}" > {output.all_Kraken_Kaiju_combined_blast_results}
            python scripts/python/reads/merge_RtPCR_kraken_kaiju_combined_report.py {output.all_Kraken_Kaiju_combined_blast_results} config/patient_metadata.tsv {params.all_Kraken_Kaiju_combined_blast_results_pcr}

            #integration of all blast files of samples
            bash scripts/shell/reads/combine_all_blast_lineage_results.sh {reads_taxa_kraken_kaiju_dir} {output.all_blast_results_virus_lineages}

            #1. Blastn_Support==Y, Percentage_Identity, Alignment_Length
            #2. virus_summary.tsv
            python scripts/python/reads/blast_confirmed_summary_generation.py {params.all_Kraken_Kaiju_combined_blast_results_pcr} {output.all_blast_results_virus_lineages} {output.virus_summary}

        """
