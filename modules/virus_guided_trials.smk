        # seq_r1 = join(preprocessing_dir, "01.trimmed/{sample}.trimmed_1.fq"),
        # seq_r2 = join(preprocessing_dir, "01.trimmed/{sample}.trimmed_2.fq"),

        # r1 = join(preprocessing_dir, "02.rmhost/{sample}_1.rmhost.fq"),
        # r2 = join(preprocessing_dir, "02.rmhost/{sample}_2.rmhost.fq"),
        #        seq_r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
        #seq_r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz"),
#This module is for guided mapping of viruses and guided assembly
checkpoint guided_virus_genomes_alignment:
    input:
        seq_r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
        seq_r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz"),
    output:
        aligned_reads = join(guided_analysis_dir, "{sample}/virus_aligned.bam"),
    conda: "base"
    threads: 10
    params:
        reference_virus_seqs = "/home/ubuntu/myMetagenomics/Trial/guided_database/Refseq_Shenzhen_dicistro-like_virus.fa",
        index_dir = join(guided_analysis_dir, "{sample}/viral_seqs_index"),
        index = join(guided_analysis_dir, "{sample}/viral_seqs_index/viral_seqs_index"),
        bam_file = join(guided_analysis_dir, "{sample}/bowtie2.bam")
    shell:
        """
        mkdir -p {params.index_dir}
        bowtie2-build {params.reference_virus_seqs} {params.index}
        bowtie2 --very-sensitive -p {threads} -x {params.index} -1 {input.seq_r1} -2 {input.seq_r2} \
        | samtools sort -O bam -@ {threads} -m 10G -o {params.bam_file}

        samtools view -b -F 4 -o {output.aligned_reads} {params.bam_file}
            
        """

rule guided_virus_genomes_assembly:
    input:
        aligned_reads = join(guided_analysis_dir, "{sample}/virus_aligned.bam"),
    output:
        assembly_dir = directory(join(guided_analysis_dir, "{sample}/megahit")),
        fasta = join(guided_analysis_dir, "{sample}/megahit/final.contigs.fa")
    conda: "base"
    params:
        reads_dir = join(guided_analysis_dir, "{sample}")
    threads: 10
    shell:
        """
        # Converts the BAM file back to FASTQ format using `samtools fastq`
        samtools fastq -1 {params.reads_dir}/potential_ref_virus_1.fastq -2 {params.reads_dir}/potential_ref_virus_2.fastq {input.aligned_reads}
        megahit -1 {params.reads_dir}/potential_ref_virus_1.fastq -2 {params.reads_dir}/potential_ref_virus_2.fastq -f -o {output.assembly_dir} -t {threads} 
        """

checkpoint diamond:
    input:
        fasta = join(guided_analysis_dir, "{sample}/megahit/final.contigs.fa"),
    output:
        matches = join(guided_analysis_dir, "{sample}/diamond/matches.tsv")
    threads: 10
    params:
        database = "/home/ubuntu/Tools/Diamond/nr_tax_full.dmnd",
        evalue = 1e-5
    benchmark:
        join(guided_analysis_dir, "{sample}/diamond/{sample}.diamond.benchmark.txt")
    conda: "diamond"
    shell:
        """
        #check if the contigs are empty
        if [ ! -s {input.fasta} ]; then
            touch {output.matches}
        else
            diamond blastx -q {input.fasta} -d {params.database} -o {output.matches} \
            --more-sensitive --evalue {params.evalue} --max-target-seqs 1 --threads {threads} \
            -f 6 qseqid sseqid pident length mismatch gapopen qcovhsp scovhsp qstart qend sstart send evalue bitscore stitle staxids sskingdoms sscinames
        fi
        """