#没有指定 rule all 的话，Snakemake 将默认执行文件中的第一个规则。但是，这可能不一定是你期望的最终目标。通过明确指定 rule all，你可以确保 Snakemake 会在执行工作流时生成你所需的目标文件。
#它将自动检查当前目录下的 Snakefile，并根据依赖关系图来确定需要执行的规则
#没有指定时想运行rmhost,命令行为：snakemake -j 2 --use-conda 2.rmhost/test1_1.rmhost.fq 2.rmhost/test1_2.rmhost.fq
#指定all后：snakemake -j 2 --use-conda  > smk.log 2>&1
# """
#  Main Snakemake workflow for: 
# Elucidating Viral Pathogens in Pediatric Unexplained Fever Using Metatranscriptomics
# """

__maintainer__ = "Zhang Tao"
__email__ = "zhangtao@hku.hk"

import os
from os.path import join
import sys
import glob
import pandas as pd
import csv
from pathlib import Path  # 推荐使用pathlib

# Directory structure
BASE_DIR = Path(".").resolve()

DATA_DIR = BASE_DIR / "resources"
Config_dir = BASE_DIR / "config"

preprocessing_dir = "00_preprocessing"
reads_taxa_kranken2_dir = "01_reads_taxa/kraken2"
reads_taxa_kaiju_dir = "01_reads_taxa/kaiju"
reads_taxa_kraken_kaiju_dir = "01_reads_taxa/kraken2_kaiju"
assembly_dir = "02_assembly"
virus_identification_dir = "03_virus_identification"
abandunce_dir = "04_abundance"
abandunce_consensus_dir = "04_abundance_consensus"
reads_contigs_dir = "05_virus_reads_contigs_integration_analysis"
virus_annotation_dir = "05_virus_annotation"
virus_component_analysis_dir = "05_virus_component_analysis"
virus_phylogenetic_dir = "06_virus_phylogenetic"

configfile: "config/config.yaml"

if not os.path.exists("logs"):
    os.makedirs("logs")

# LOAD METADATA
metadata = pd.read_csv(Config_dir / "samples.tsv", sep="\t", index_col=0)
SAMPLES = metadata["name"].str.strip().tolist()

# Create a dictionary with sample IDs as keys and file paths as values
names = metadata["name"].str.strip().to_dict()

# Define a function to construct the file path using the sample name
def get_file_path(sample_id, name, read_number, base_dir=DATA_DIR):
     sample_name = name + "_" + read_number + ".fastq.gz"
     return os.path.join(base_dir, sample_id, sample_name)
    
all_outfiles = [
     # preprocessing
     # #使用expand来自动组合出所有的目标文件
     #expand(join(preprocessing_dir, "00.raw_qc/fastqc/{sample}_{R}_fastqc.html"), sample=names.keys(), R=["1","2"]),
     # expand(join(preprocessing_dir, "00.raw_qc/trim_fastqc/{sample}_{R}_fastqc.html"), sample=names.keys(), R=["1","2"]),
     # # join(reads_taxa_kranken2_dir, "PMT8346_trim", "PMT8346.b.krona.txt"),
     # expand(join(reads_taxa_kranken2_dir, "{sample}", "{sample}.b.krona.txt"), sample=names.keys()),
     # expand(join(reads_taxa_kranken2_dir, "{sample}", "{sample}.b.krona.txt"), sample=names.keys()),
     # directory(join(reads_taxa_kraken_kaiju_dir, "PMT5836/blast_results")),
     # join(virus_identification_dir,"all_combined_virus_results_lineages.tsv"),
     # join(reads_taxa_kraken_kaiju_dir, "virus_summary.tsv"),
     # join(reads_contigs_dir,"virus_summary_with_contigs.tsv")
     join(reads_contigs_dir,"virus_summary_with_contigs_indexhopping.filtered.final.dnarna.tsv")
     #join(abandunce_dir, "all_concatenated_merged_abundance_blast_taxa.no.phages_with_genome.txt")
     # expand(join(abandunce_dir, "{sample}/virus_abundance.txt"), sample=names.keys())

     #expand(join(virus_identification_dir, "{sample}/blastn/matches.tsv"), sample=names.keys()),
     #expand(join(virus_identification_dir, "{sample}", "diamond_blastn/combined_virus_results.txt"), sample=names.keys()),
     #expand(join(virus_identification_dir, "{sample}/blastn/potential_non_viruses.tsv"), sample=names.keys()),
     #expand(join(virus_identification_dir, "{sample}/blastn/potential_endogenous_viruses.tsv"), sample=names.keys()),
     #expand(join(abandunce_dir, "{sample}", "merged_abundance_blast_taxa.txt"), sample=names.keys()),
     # join(virus_identification_dir, "all_combined_virus_results.no.phages.txt")
     # expand(join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq"), sample=names.keys()),

     # join(abandunce_dir, "quantification_result.txt"),
     # join(virus_component_analysis_dir, "quantification_rpm_blast_results.tsv"),
     # join(virus_component_analysis_dir, "quantification_result.txt.no.phages.RPM_length_filtered")
     # join(virus_component_analysis_dir, "OK.txt")
]

rule all:
   input: all_outfiles

# test var
# var = get_file_path("kkk","name","2")
# rule print_var:
#     run:
#         print(var)
# rule print_var_shell:
#     shell:
#         """
#         echo '{var}'
#         """
# include: "modules_pipeline/preprocessing.smk"
# include: "modules_pipeline/reads_taxa_assignment.smk"
# include: "modules_pipeline/assembly.smk"
# include: "modules_pipeline/virus_identification.smk"

# include: "modules_pipeline/virus_abundance_calculation.smk"
# include: "modules_pipeline/novel_virus_identification.smk"
# include: "modules_pipeline/virus_annotation.smk"
# include: "modules_pipeline/Bacterial_fungal_identification.smk"
# include: "modules_pipeline/phylogenetic_analysis.smk"

include: "modules/preprocessing.smk"
#include: "modules/reads_taxa_assignment.smk"
include: "modules/assembly.smk"
include: "modules/virus_identification_reads_based.smk"
include: "modules/virus_identification_contig_based.smk"
include: "modules/virus_abundance_calculation.smk"
include: "modules/virus_filter.smk"

#include: "modules/virus_component_analysis.smk"
# include: "modules/virus_phylogenetic.smk"



