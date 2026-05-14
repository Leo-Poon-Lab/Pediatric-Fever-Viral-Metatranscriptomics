checkpoint megahit:
    input:
        r1 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_1.fq.gz"),
        r2 = join(preprocessing_dir, "03.rmRrna/{sample}_rmhost_rmRrna_2.fq.gz")
    output:
        output_dir = directory(join(assembly_dir, "{sample}/megahit")),
        fasta = join(assembly_dir, "{sample}/megahit/final.contigs.fa")
    threads: 10
    benchmark:
        join(assembly_dir, "{sample}/megahit/{sample}.megahit.benchmark.txt")
    conda: "../envs/workflow.yaml"
    shell:
        "megahit -1 {input.r1} -2 {input.r2} -f -o {output.output_dir} -t {threads}"