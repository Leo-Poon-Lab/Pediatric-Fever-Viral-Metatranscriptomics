checkpoint diamond:
    input:
        contigs = join(assembly_dir, "{sample}/megahit/final.contigs.fa"),
    output:
        matches = join(virus_identification_dir, "{sample}/diamond/matches.tsv")
    threads: 10
    params:
        database = Path(config["databases"]["protein"]["nr"]),
        evalue = 1e-5,
        matches_header = join(virus_identification_dir, "{sample}/diamond/matches_header.tsv")
    benchmark:
        join(virus_identification_dir, "{sample}/diamond/{sample}.diamond.benchmark.txt")
    conda: "diamond"
    shell:
        """
        check if the contigs are empty
        if [ ! -s {input.contigs} ]; then
            touch {output.matches}
        else
            diamond blastx -q {input.contigs} -d {params.database} -o {output.matches} \
            --more-sensitive --evalue {params.evalue} --max-target-seqs 1 --threads {threads} \
            -f 6 qseqid sseqid pident length mismatch gapopen qcovhsp scovhsp qstart qend sstart send evalue bitscore stitle staxids sskingdoms sscinames

            # Define header
            echo -e "Query_Sequence_ID\tSubject_Sequence_ID\tPercentage_of_Identical_Matches\tAlignment_Length\tNumber_of_Mismatches\tNumber_of_Gap_Openings\tQuery_Coverage_Per_HSP\tSubject_Coverage_Per_HSP\tStart_of_Alignment_in_Query\tEnd_of_Alignment_in_Query\tStart_of_Alignment_in_Subject\tEnd_of_Alignment_in_Subject\tE_Value\tBit_Score\tSubject_Title\tSubject_Taxonomy_IDs\tSubject_Kingdoms\tSubject_Scientific_Name" >  {params.matches_header}

            cat {params.matches_header} >> {params.matches_header}

        fi
        """

rule identify_potential_virus:
    input:
        matches = join(virus_identification_dir, "{sample}/diamond/matches.tsv"),
        contigs = join(assembly_dir, "{sample}/megahit/final.contigs.fa"),
    output:
        virus_matches = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs_all.tsv"),
        virus_matches_contigs = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs.fa"),
        virus_no_phage_matches = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs_all_no_phage.tsv"),
    params:
        matches_dir = join(virus_identification_dir, "{sample}/diamond")
    conda: "base"
    shell:
        """
        #check if the matches are empty
        if [ ! -s {input.matches} ]; then
            touch {output.virus_matches}
            touch {output.virus_no_phage_matches}
            touch {output.virus_matches_contigs}
        else
            # generate formal taxonomic table
            bash scripts/shell/contigs/taxonkit_retrive_taxa_from_taxaid.sh {params.matches_dir} {params.matches_dir}

            # Ensure the output file from taxonkit exists before proceeding
            if [ -s {output.virus_matches} ]; then

                if grep -q -Evi "phage|bacteriophage|Caudovir|Microviridae|Microvirus|Myoviridae|Prokaryotic|Siphoviridae|Siphovirus|Podoviridae|Podovirus|Inoviridae|leviviridae|Herelleviridae" {output.virus_matches}; then
                # Non-phage content found
                    grep -Evi "phage|bacteriophage|Caudovir|Microviridae|Microvirus|Myoviridae|Prokaryotic|Siphoviridae|Siphovirus|Podoviridae|Podovirus|Inoviridae|leviviridae|Herelleviridae" {output.virus_matches} > {output.virus_no_phage_matches}
                    seqkit grep -f <(cut -f1 {output.virus_no_phage_matches}) {input.contigs} > {output.virus_matches_contigs}
                else
                    touch {output.virus_no_phage_matches}
                    touch {output.virus_matches_contigs}
                fi
                
            else
                touch {output.virus_matches}
                touch {output.virus_no_phage_matches}
                touch {output.virus_matches_contigs}
            fi           
        fi
        
        """

checkpoint blastn:
    input:
        virus_matches_contigs = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs.fa")
    output:
        matches = join(virus_identification_dir, "{sample}/blastn/matches.tsv")
    threads: 10
    params:
        database = Path(config["databases"]["nucleotide"]["nt"]),
        task = "blastn",
        evalue = 1e-10,
        matches_dir = join(virus_identification_dir, "{sample}/blastn"),
        matches_header = join(virus_identification_dir, "{sample}/blastn/matches_header.tsv")
    benchmark:
        join(virus_identification_dir, "{sample}/blastn/{sample}.blastn.benchmark.txt")
    conda: "base"
    shell:
        """
        #check if the contigs are empty
        if [ ! -s {input.virus_matches_contigs} ]; then
            touch {output.matches}
        else
            blastn -query {input.virus_matches_contigs} -db {params.database} -task {params.task} -outfmt "\
            6 qseqid qlen sseqid slen pident nident ppos positive length qcovs qcovhsp \
            qcovus mismatch gapopen qstart qend sstart send gaps evalue bitscore sstrand \
            staxid stitle" -max_target_seqs 3 \
            -evalue {params.evalue} -num_threads {threads} -out {output.matches}

            # Add header
            echo -e "Query_ID\tQuery_Length\tSubject_ID\tSubject_Length\tPercentage_Identity\tNum_Identical_Matches\tPercentage_Positive\tNum_Positive_Matches\tAlignment_Length\tQuery_Coverage_Percent\tQuery_Coverage_Percent_HSP\tQuery_Coverage_Percent_Unified\tNum_Mismatches\tNum_Gap_Openings\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tNum_Gaps\tE_value\tBit_Score\tSubject_Strand\tSubject_Taxonomy_ID\tSubject_Title" > {params.matches_header}
            cat {output.matches} >> {params.matches_header}
        
        fi
        """

rule identify_non_virus_blastn:
    input:
        virus_matches_blastn = join(virus_identification_dir, "{sample}/blastn/matches.tsv"),
    output:
        non_virus = join(virus_identification_dir, "{sample}/blastn/potential_non_viruses.tsv"),
        endogenous_virus = join(virus_identification_dir, "{sample}/blastn/potential_endogenous_viruses.tsv")
    params:
        blastn_dir = join(virus_identification_dir, "{sample}/blastn"),
    conda: "base"
    shell:
        """
        if [ ! -s {input.virus_matches_blastn} ]; then
            touch {output.non_virus}
            touch {output.endogenous_virus}
        else
            bash scripts/shell/contigs/taxonkit_retrive_taxa_from_taxaid_blastn.sh {params.blastn_dir} {output.non_virus} {output.endogenous_virus}
        fi
        """

# combine blastx and blastn results
rule anlysis_blastx_blastn_virus:
    input:
        virus_matches_diamond = join(virus_identification_dir, "{sample}/diamond/matches.tsv"),
        virus_matches_blastn = join(virus_identification_dir, "{sample}/blastn/matches.tsv"),
    output:
        combined_virus_results = join(virus_identification_dir, "{sample}/diamond_blastn/combined_virus_results.txt")
    conda: "base"
    shell:
        """
        if [ ! -s {input.virus_matches_diamond} ]; then
            touch {output.combined_virus_results}
        else
            python scripts/python/contigs/combine_blast_results.py {virus_identification_dir} {wildcards.sample} {input.virus_matches_diamond} {input.virus_matches_blastn} {output.combined_virus_results}
        fi
        
        """
        
# Waiting for all samples finished
# integrate combined blastx and blastn results for all samples
rule combine_blastx_blastn_virus:
    input:
        work_dir = virus_identification_dir,
        # 明确列出所有样本的中间文件
        combined_per_sample = expand(
            join(virus_identification_dir, "{sample}/diamond_blastn/combined_virus_results.txt"),
            sample=SAMPLES
        )
    output:
        all_combined_diamond_blastn = join(virus_identification_dir, "all_combined_virus_results.txt"),
        all_combined_diamond_blastn_lineage = join(virus_identification_dir, "all_combined_virus_results_lineages.tsv"),
        all_combined_diamond_blastn_without_phage = join(virus_identification_dir, "all_combined_virus_results.no.phages.txt"),
        all_combined_diamond_blastn_without_phage_lineage = join(virus_identification_dir, "all_combined_virus_results.no.phages_lineages.tsv")
    conda: "base"
    shell:
        """ 
        bash scripts/shell/contigs/combine_all_diamond_blastn_results.sh {input.work_dir} {output.all_combined_diamond_blastn}

        grep -Evi "phage|bacteriophage|Caudovir|Microviridae|Microvirus|Myoviridae|Prokaryotic|Siphoviridae|Siphovirus|Podoviridae|Podovirus|Inoviridae|leviviridae|Herelleviridae|Ackermannviridae|Crassvirales"  {output.all_combined_diamond_blastn} > {output.all_combined_diamond_blastn_without_phage}

        #Add Lineages info
        bash scripts/shell/contigs/add_lineages.sh {output.all_combined_diamond_blastn} {output.all_combined_diamond_blastn_lineage} {output.all_combined_diamond_blastn_without_phage} {output.all_combined_diamond_blastn_without_phage_lineage}

        if [ ! -s {output.all_combined_diamond_blastn} ]; then
            touch {output.all_combined_diamond_blastn}
            touch {output.all_combined_diamond_blastn_lineage}
        fi
        """