configfile: "config/config.yaml"


SAMPLES = config["samples"]


rule all:
    input:
        "results/differential_expression/DE_results.csv",
        "results/multiqc/multiqc_report.html",
        expand(
            "results/alignments/{sample}.sorted.bam.bai",
            sample=SAMPLES
        )

# -------------------------
# Fast QC (Pre-Trimming)
# -------------------------
rule fastqc_raw:
    input:
        r1="data/raw/{sample}_1.fastq.gz",
	r2="data/raw/{sample}_2.fastq.gz"
    output:
        html1="results/qc/raw/{sample}_1_fastqc.html",
        zip1="results/qc/raw/{sample}_1_fastqc.zip",
        html2="results/qc/raw/{sample}_2_fastqc.html",
        zip2="results/qc/raw/{sample}_2_fastqc.zip"
    threads: 4
    shell:
        """
        mkdir -p results/qc/raw
        fastqc {input.r1} {input.r2} \
        -o results/qc/raw \
        --threads {threads}
        """


# -------------------------
# Trimming 
# -------------------------
rule trim:
    input:
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz"
    output:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz",
        unp1="results/trimmed/{sample}_R1.unpaired.fastq.gz",
        unp2="results/trimmed/{sample}_R2.unpaired.fastq.gz"
    threads: 8
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
	    -phred33 \
            {input.r1} {input.r2} \
            {output.r1} {output.unp1} \
            {output.r2} {output.unp2} \
            SLIDINGWINDOW:4:20 MINLEN:50
        """
# -------------------------
# Fast QC (Post-Trimming)
# -------------------------
rule fastqc_trimmed:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        html1="results/qc/trimmed/{sample}_R1.trimmed_fastqc.html",
        html2="results/qc/trimmed/{sample}_R2.trimmed_fastqc.html",
        zip1="results/qc/trimmed/{sample}_R1.trimmed_fastqc.zip",
        zip2="results/qc/trimmed/{sample}_R2.trimmed_fastqc.zip"
    threads: 4
    shell:
        """
        mkdir -p results/qc/trimmed

        fastqc {input.r1} {input.r2} \
            --outdir results/qc/trimmed \
            --threads {threads}
        """

# -------------------------
# HISAT2 index
# -------------------------
rule hisat2_index:
    input:
        config["reference"]["genome"]
    output:
        directory("data/reference/hisat2_index")
    threads: 4
    log:
        "logs/hisat2_index.log"
    shell:
        """
        mkdir -p data/reference/hisat2_index logs

        hisat2-build {input} data/reference/hisat2_index/mouse_index \
        > {log} 2>&1
        """


# -------------------------
# HISAT2 alignment
# -------------------------
rule align:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz", 
	index=rules.hisat2_index.output
    output:
        bam="results/alignments/{sample}.sorted.bam",
	bai="results/alignments/{sample}.sorted.bam.bai"
    threads: 4
    shell:
        """
        hisat2 \
            -x data/reference/hisat2_index/mouse_index \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} \
        | samtools sort \
            -@ {threads} \
            -o {output.bam}

	samtools index {output.bam}
        """

# -------------------------
# flagstat alignment QC
# -------------------------

rule flagstat:
    input:
        "results/alignments/{sample}.sorted.bam"
    output:
        "results/qc/alignment/{sample}.flagstat.txt"
    log:
        "logs/flagstat/{sample}.log"
    shell:
        """
        mkdir -p results/qc/alignment logs/flagstat

        samtools flagstat {input} > {output}

        echo "Done" > {log}
        """

# -------------------------
# featureCounts (fixed BAM input handling)
# -------------------------
rule featurecounts:
    input:
        bams=expand("results/alignments/{sample}.sorted.bam", sample=SAMPLES),
        gtf=config["reference"]["gtf"]
    output:
        "results/counts/gene_counts.txt",
        "results/counts/gene_counts.txt.summary"
    threads: 4
    shell:
        """
        mkdir -p results/counts logs

        featureCounts -T {threads} \
	-p \
        -a {input.gtf} \
        -o results/counts/gene_counts.txt \
        -t exon \
        -g gene_id \
        {input.bams}
        """


# -------------------------
# deseq2 Analysis
# -------------------------
rule deseq2:
    input:
        counts="results/counts/gene_counts.txt",
	metadata="metadata/samples.csv"
    output:
        "results/differential_expression/DE_results.csv"
    log:
        "logs/deseq2.log"
    script:
        "scripts/deseq2_analysis.R"

# -------------------------
# MultiQC
# -------------------------
rule multiqc:
    input:
        expand("results/qc/raw/{sample}_1_fastqc.zip", sample=SAMPLES),
        expand("results/qc/raw/{sample}_2_fastqc.zip", sample=SAMPLES),
        expand("results/qc/trimmed/{sample}_R1.trimmed_fastqc.zip", sample=SAMPLES),
        expand("results/qc/trimmed/{sample}_R2.trimmed_fastqc.zip", sample=SAMPLES),
        expand("results/qc/alignment/{sample}.flagstat.txt", sample=SAMPLES),
        "results/counts/gene_counts.txt"
    output:
        "results/multiqc/multiqc_report.html"
    shell:
        """
        mkdir -p results/multiqc logs

        multiqc results \
            -o results/multiqc \
            > logs/multiqc.log 2>&1
        """
