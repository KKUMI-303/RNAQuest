configfile: "config/config.yaml"


rule all:
    input:
        "results/differential_expression/DE_results.csv"


rule fastqc:
    input:
        "data/raw/{sample}.fastq"
    output:
        "results/qc/{sample}_fastqc.html"
    threads: 2
    log:
        "logs/fastqc/{sample}.log"
    shell:
        """
        fastqc {input} -o results/qc/ > {log} 2>&1
        """

rule trim:
    input:
        r1="data/raw/{sample}_R1.fastq",
        r2="data/raw/{sample}_R2.fastq"
    output:
        r1="data/trimmed/{sample}_R1_paired.fastq",
        r2="data/trimmed/{sample}_R2_paired.fastq",
        r1_unpaired="data/trimmed/{sample}_R1_unpaired.fastq",
        r2_unpaired="data/trimmed/{sample}_R2_unpaired.fastq"
    threads: 4
    log:
        "logs/trimmomatic/{sample}.log"
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} \
        {output.r1} {output.r1_unpaired} \
        {output.r2} {output.r2_unpaired} \
        SLIDINGWINDOW:4:20 MINLEN:50 \
        > {log} 2>&1
        """


rule star_index:
    input:
        genome=config["reference"]["genome"],
        gtf=config["reference"]["gtf"]
    output:
        directory("data/reference/star_index")
    threads: 4
    log:
        "logs/star_index.log"
    shell:
        """
        STAR --runMode genomeGenerate \
        --genomeDir data/reference/star_index \
        --genomeFastaFiles {input.genome} \
        --sjdbGTFfile {input.gtf} \
        --runThreadN {threads} \
        > {log} 2>&1
        """



rule align:
    input:
        r1="data/trimmed/{sample}_R1_paired.fastq",
        r2="data/trimmed/{sample}_R2_paired.fastq",
        index="data/reference/star_index"
    output:
        "results/alignments/{sample}_Aligned.out.sam"
    threads: 4
    log:
        "logs/star/{sample}.log"
    shell:
        """
        STAR --genomeDir {input.index} \
        --readFilesIn {input.r1} {input.r2} \
        --outFileNamePrefix results/alignments/{wildcards.sample}_ \
        --outSAMtype SAM \
        --runThreadN {threads} \
        > {log} 2>&1
        """



rule bam_sort:
    input:
        "results/alignments/{sample}_Aligned.out.sam"
    output:
        bam="results/alignments/{sample}.sorted.bam"
    threads: 2
    log:
        "logs/samtools/{sample}.log"
    shell:
        """
        samtools view -Sb {input} | samtools sort -o {output.bam}
        samtools index {output.bam} > {log} 2>&1
        """

rule featurecounts:
    input:
        bam=expand(
            "results/alignments/{sample}.sorted.bam",
            sample=config["samples"]
        ),
        gtf=config["reference"]["gtf"]
    output:
        "results/counts/gene_counts.txt"
    threads: 4
    log:
        "logs/featurecounts.log"
    shell:
        """
        featureCounts -T {threads} -a {input.gtf} -o {output} {input.bam} \
        > {log} 2>&1
        """

rule deseq2:
    input:
        "results/counts/gene_counts.txt"
    output:
        "results/differential_expression/DE_results.csv"
    log:
        "logs/deseq2.log"
    script:
        "../scripts/deseq2_analysis.R"
