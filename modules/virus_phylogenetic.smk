##################################################################
# 选择病毒进行系统发育分析
# 看大框架树和看精细树
# 看不同病毒的genome coverage情况：all_mapping_summary.xlsx
# 看不同病毒的Identity情况，判断是否是新病毒：virus_summary_with_contigs_indexhopping.final.filtered.added.dna-rna.host.updated.animal.merged
##################################################################

rule mafft:
    input:
        seqs = join(virus_phylogenetic_dir, "Dicistro-like_virus/samples_replication_protein_sequences.fasta"),
    output:
        filtered_contigs = join(virus_phylogenetic_dir, "Dicistro-like_virus/samples_replication_protein_sequences.maffs.fas"),
    threads: 12
    params:
        length = 300
    conda: "phylogeny"
    shell:
        """ 
        #先ref，再add
        mafft --auto Dicistroviridae_reference.fa > reference_alignment.maff.fas

        mafft --auto  --thread 12 --adjustdirection --add query_protein_sequences.fasta --reorder reference_alignment.maff.fas > alignment.maffs.fas

        mafft --auto  --thread 12 --adjustdirection --reorder samples_replication_protein_sequences.fasta > samples_replication_protein_sequences.maffs.fas
        """

rule trimAI:
    input:
        seqs = join(virus_phylogenetic_dir, "Dicistro-like_virus/alignment.maffs.fas"),
    output:
        filtered_contigs = join(virus_phylogenetic_dir, "Dicistro-like_virus/alignment.maffs.trimal.fas"),
    threads: 12
    params:
        length = 300
    conda: "phylogeny"
    shell:
        """ 
        trimal -in alignment.maffs.fas -out alignment.maffs.trimal.fas -automated1
        """

rule iqtree:
    input:
        seqs = join(virus_phylogenetic_dir, "Dicistro-like_virus/alignment.maffs.trimal.fas"),
    output:
        filtered_contigs = join(virus_phylogenetic_dir, "Dicistro-like_virus/alignment.maffs.trimal.iqtree"),
    threads: 12
    params:
        length = 300
    conda: "phylogeny"
    shell:
        """ 
        mkdir iqtree
        mkdir log
        
        iqtree2 -T 10 -s alignment.maffs.trimal.fas --prefix ./iqtree/Dicistro-like -m MFP -B 1000 -bnni -alrt 1000 >log/iqtree.log 2>&1

        """

rule visiualize:
    input:
        seqs = join(virus_phylogenetic_dir, "Dicistro-like_virus/alignment.maffs.trimal.iqtree"),
    output:
        filtered_contigs = join(virus_phylogenetic_dir, "Dicistro-like_virus/visiualization/alignment.maffs.trimal.iqtree.treefile"),
    threads: 12
    params:
        length = 300
    conda: "phylogeny"
    shell:
        """ 
        # visiualize
        iTOL -s ./iqtree/Dicistro-like.treefile -t ./iqtree/Dicistro-like.treefile -o ./iqtree/Dicistro-like_iTOL -w 1000 -x 1000 -y 1000 -z 1000

        # https://link.springer.com/article/10.1186/s13073-025-01447-3#Sec16
        # 两种tree
        Rscript /home/ubuntu/myMetagenomics/Trial/06_virus_phylogenetic/Enterovirus_A/global_tree_visualiization.R
        Rscript /home/ubuntu/myMetagenomics/Trial/06_virus_phylogenetic/Dicistro-like_virus/visiualization/global_tree_visualiization.R
        """


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
#                 touch {output.virus_matches_contigs}
#             fi
#         fi
        
#         """
