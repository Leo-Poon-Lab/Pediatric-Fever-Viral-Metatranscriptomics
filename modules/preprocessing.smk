rule raw_fastqc:
    input:
        r1 = lambda wildcards: get_file_path(wildcards.sample,names[wildcards.sample],"1"),
        r2 = lambda wildcards: get_file_path(wildcards.sample,names[wildcards.sample],"2")
    output:
        r1_qc = join(preprocessing_dir, "00.raw_qc/fastqc/{sample}_1_fastqc.html"),
        r2_qc = join(preprocessing_dir, "00.raw_qc/fastqc/{sample}_2_fastqc.html")
    params:
        outdir = directory(join(preprocessing_dir, "00.raw_qc/fastqc"))
    threads: 10
    conda:
        "envs/workflow.yaml"
    shell:
        """
        fastqc {input.r1} -t {threads} --outdir {params.outdir} - > {output.r1_qc}
        fastqc {input.r2} -t {threads} --outdir {params.outdir} - > {output.r2_qc}
        """

rule trim_fastp:
    input:
        r1 = lambda wildcards: get_file_path(wildcards.sample,names[wildcards.sample],"1"),
        r2 = lambda wildcards: get_file_path(wildcards.sample,names[wildcards.sample],"2")
    output:
        trim_r1=join(preprocessing_dir, "01.trimmed/{sample}.trimmed_1.fq"),
        trim_r2=join(preprocessing_dir, "01.trimmed/{sample}.trimmed_2.fq"),
    threads: 10
    params:
        min_len = 36,
        qualified_quality_phred = 20,
    log:
        json = join(preprocessing_dir, "01.trimmed/{sample}.fastp.json"),
        html = join(preprocessing_dir, "01.trimmed/{sample}.fastp.html"),
        fastp_log = join(preprocessing_dir, "01.trimmed/{sample}.fastp.log")
    conda:
        "envs/workflow.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.trim_r1} -O {output.trim_r2} \
        -5 -3 -c --detect_adapter_for_pe \
        --qualified_quality_phred {params.qualified_quality_phred} -w {threads} --length_required {params.min_len} \
        -j {log.json} -h {log.html} 2> {log.fastp_log}
        """

rule trim_qc:
    input:
        trim_r1=join(preprocessing_dir, "01.trimmed/{sample}.trimmed_1.fq"),
        trim_r2=join(preprocessing_dir, "01.trimmed/{sample}.trimmed_2.fq"),
    output:
        trim_r1_qc = join(preprocessing_dir, "00.raw_qc/trim_fastqc/{sample}_1_fastqc.html"),
        trim_r2_qc = join(preprocessing_dir, "00.raw_qc/trim_fastqc/{sample}_2_fastqc.html")
    params:
        outdir = directory(join(preprocessing_dir, "00.raw_qc/trim_fastqc"))
    threads: 10
    conda:
        "envs/workflow.yaml"
    shell:
        """
        fastqc {input.trim_r1} -t {threads} --outdir {params.outdir} - > {output.trim_r1_qc}
        fastqc {input.trim_r2} -t {threads} --outdir {params.outdir} - > {output.trim_r2_qc}
        """

rule rmhost:
    input:
        trim_r1=join(preprocessing_dir, "01.trimmed/{sample}.trimmed_1.fq"),
        trim_r2=join(preprocessing_dir, "01.trimmed/{sample}.trimmed_2.fq"),
    output:
        rmhost_r1 = join(preprocessing_dir, "02.rmhost/{sample}_1.rmhost.fq"),
        rmhost_r2 = join(preprocessing_dir, "02.rmhost/{sample}_2.rmhost.fq")
    params:
        index_bwa2 = Path(config["databases"]["human"]["bwa_index"]),
        aligned_reads = join(preprocessing_dir, "02.rmhost/{sample}.aligned.sam"),
    threads: 10
    conda:
        "envs/workflow.yaml"
    shell:
        """
        bwa-mem2 mem -t {threads} -M {params.index_bwa2} {input.trim_r1} {input.trim_r2} -o {params.aligned_reads}
        samtools fastq -f 12 -F 256 {params.aligned_reads} -1 {output.rmhost_r1} -2 {output.rmhost_r2}
        rm {params.aligned_reads}
        """

rule rmRrna:
    input:
        rmhost_r1=join(preprocessing_dir, "02.rmhost/{sample}_1.rmhost.fq"),
        rmhost_r2=join(preprocessing_dir, "02.rmhost/{sample}_2.rmhost.fq"),

    output:
        rmRrna_r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
        rmRrna_r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz")
    params:
        log = join(preprocessing_dir, "03.rmRrna/{sample}_rRNA.log"),
        workdir = Path(config["databases"]["rrna"]["workdir"]),
        outputRrna =  join(preprocessing_dir, "03.rmRrna/{sample}_rRNA"),
        outputnon_Rrna = join(preprocessing_dir, "03.rmRrna/{sample}_non_rRNA"),

        db_archaea = Path(config["databases"]["rrna"]["sortmerna"]["archaea_5s"]),
        db_bacteria = Path(config["databases"]["rrna"]["sortmerna"]["bacteria_5s"]),
        db_eukaryota = Path(config["databases"]["rrna"]["sortmerna"]["eukaryota_5s"]),
        db_lsu = Path(config["databases"]["rrna"]["sortmerna"]["lsu"]),
        db_ssu = Path(config["databases"]["rrna"]["sortmerna"]["ssu"])
    threads: 20
    conda:
        "base"
    shell:
        """
        if [ "$(ls -A {params.workdir}/kvdb)" ]; then
            echo "Files found in kvdb directory. Removing them..."
            rm {params.workdir}/kvdb/* 
        fi
        sortmerna --index 0 --threads {threads} \
        --ref {params.db_archaea} \
        --ref {params.db_bacteria} \
        --ref {params.db_eukaryota} \
        --ref {params.db_lsu} \
        --ref {params.db_ssu} \
        --workdir {params.workdir} \
        --reads {input.rmhost_r1} \
        --reads {input.rmhost_r2} \
        --aligned {params.outputRrna} --other {params.outputnon_Rrna} \
        --paired_in -out2 --fastx  > {params.log} 2>&1

        pigz -p {threads} {params.outputnon_Rrna}_fwd.fq
        pigz -p {threads} {params.outputnon_Rrna}_rev.fq

        mv {params.outputnon_Rrna}_fwd.fq.gz {output.rmRrna_r1}
        mv {params.outputnon_Rrna}_rev.fq.gz {output.rmRrna_r2}
    
        """

#rule rmRrna_bowtie2:
#     input:
#         rmhost_r1=join(preprocessing_dir, "02.rmhost/{sample}_1.rmhost.fq"),
#         rmhost_r2=join(preprocessing_dir, "02.rmhost/{sample}_2.rmhost.fq"),
#        
#     output:
#         rmRrna_r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
#         rmRrna_r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz")
#        
#     params:
#         index = Path(config["databases"]["rrna"]["bowtie2"]["rRNA_db"]),
#         outputRrna =  join(preprocessing_dir, "03.rmRrna/{sample}_rRNA"),
#         outputnon_Rrna = join(preprocessing_dir, "03.rmRrna/{sample}_non_rRNA")
#     threads: 20
#     conda:
#         "base"
#     shell:
#         """
#         bowtie2 -p {threads} -x {params.index} -1 {input.rmhost_r1} -2 {input.rmhost_r2} \
#         --un-conc {params.outputnon_Rrna}.fq -S /dev/null > {params.outputRrna}.log 2>&1
#
#         mv {params.outputnon_Rrna}.1.fq {output.rmRrna_r1}
#         mv {params.outputnon_Rrna}.2.fq {output.rmRrna_r2}
#    
#         """