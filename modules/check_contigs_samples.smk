rule assembly_quast:
    input:
        contigs = join(assembly_dir, "{sample}/megahit/final.contigs.fa"),
    output:
        output_dir = directory(join(assembly_dir, "{sample}/megahit/quast_results")),
    params:
        threads = 5
    conda: "base"
    shell:
        "python /home/ubuntu/anaconda3/bin/quast.py {input.contigs} --min-contig 0 -t {params.threads} -o {output.output_dir} --silent --no-plots"

# rule create_empty_contigs:
#     output:
#         contigs = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs.fa")
#     shell:
#         """
#         # 全部覆盖了
#                 touch {output.contigs}
#         """
# select kingdom virus as potential virus
rule identify_potential_virus_re:
    input:
        matches = join(virus_identification_dir, "{sample}/diamond/matches.tsv"),
        contigs = join(assembly_dir, "{sample}/megahit/final.contigs.fa"),
    output:
        virus_matches = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs_all.tsv"),
        virus_matches_contigs = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs.fa")
    params:
        matches_dir = join(virus_identification_dir, "{sample}/diamond")
    conda: "base"
    shell:
        """
        #check if the matches are empty
        if [ ! -s {input.matches} ]; then
            touch {output.virus_matches}
            touch {output.virus_matches_contigs}
        else
            # generate formal taxonomic table
            # bash /home/ubuntu/Tools/Scripts/taxonkit_retrive_taxa_from_taxaid.sh {params.matches_dir}
            bash /home/ubuntu/Tools/Scripts/taxonkit_retrive_taxa_from_taxaid_try.sh {params.matches_dir}  
            if grep -qi "Virus" {output.virus_matches}; then
                seqkit grep -f <(cut -f1 {output.virus_matches}) {input.contigs} > {output.virus_matches_contigs}
            else
                touch {output.virus_matches}
                touch {output.virus_matches_contigs}
            fi
        fi
        echo "fininished"
        """
#由于改动了taxonkit_retrive_taxa_from_taxaid.sh，所以需要重新生成potential_virus_contigs_all.tsv，看看是否有变化
#单纯检验而已，如果有新的测序数据，不需要跑这个rule
rule identify_potential_virus_re_2:
    input:
        matches = join(virus_identification_dir, "{sample}/diamond/matches.tsv"),
        contigs = join(assembly_dir, "{sample}/megahit/final.contigs.fa"),
    output:
        virus_matches = join(virus_identification_dir, "{sample}/diamond/try/potential_virus_contigs_all.tsv"),
        virus_no_phage_matches = join(virus_identification_dir, "{sample}/diamond/try/potential_virus_contigs_all_no_phage.tsv"),
        virus_matches_contigs = join(virus_identification_dir, "{sample}/diamond/try/potential_virus_contigs.fa")
    params:
        old_matches = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs_all.tsv"),
        matches_dir = join(virus_identification_dir, "{sample}/diamond"),
        output_dir = join(virus_identification_dir, "{sample}/diamond/try")
    conda: "base"
    shell:
        """
        #check if the matches are empty
        if [ ! -s {input.matches} ]; then
            touch {output.virus_matches}
            touch {output.virus_matches_contigs}
            touch {output.virus_no_phage_matches}
        else

            bash /home/ubuntu/Tools/Scripts/taxonkit_retrive_taxa_from_taxaid_try.sh {params.matches_dir} {params.output_dir}
            
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
                
                # Compare old and new virus matches files if the old file exists
                if [ -s {params.old_matches} ]; then
                    if diff "{params.old_matches}" "{output.virus_no_phage_matches}"; then
                        echo "No differences found."
                    else
                        echo "Differences found between {params.old_matches} and {output.virus_no_phage_matches}." >> /home/ubuntu/myMetagenomics/Trial/03_virus_identification/diff_virus_matches.txt
                    fi
                fi
            else
                touch {output.virus_matches}
                touch {output.virus_no_phage_matches}
                touch {output.virus_matches_contigs}
            fi                
        fi
        """

rule virus_quast:
    input:
        contigs = join(virus_identification_dir, "{sample}/diamond/potential_virus_contigs.fa"),
    output:
        output_dir = directory(join(virus_identification_dir, "{sample}/diamond/quast_results")),
    params:
        threads = 5,
        matches_dir = join(virus_identification_dir, "{sample}/diamond"),
        try_dir = join(virus_identification_dir, "{sample}/diamond/try")
    conda: "base"
    shell:
        """
        cp {params.try_dir}/* {params.matches_dir}
        
        if [ ! -s {input.contigs} ]; then
            mkdir -p {output.output_dir}
        else
            python /home/ubuntu/anaconda3/bin/quast.py {input.contigs} --min-contig 0 -t {params.threads} -o {output.output_dir} --silent --no-plots
        fi
        
        """
