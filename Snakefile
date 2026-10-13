configfile: "config/config.yaml"
RAW = config["paths"]["raw"]
RESULTS = config["paths"]["results"]
LOGS = config["paths"]["logs"]
HISAT2_INDEX = config["reference"]["hisat2_index"]
META = config["paths"]["metadata"]
SCRIPTS = config["scripts"]["test_script"]
import csv

with open(META) as f:
    reader = csv.DictReader(f)
    SAMPLES = [row["sample"] for row in reader]
rule all:
    input:
        expand(
            f"{RESULTS}/alignments/{{sample}}.sorted.bam.bai",
            sample=SAMPLES
        ),
	f"{RESULTS}/counts/gene_counts.txt",
        f"{RESULTS}/multiqc/multiqc_report.html",
        f"{RESULTS}/differential_expression/DE_results.csv",
        f"{RESULTS}/differential_expression/comparison.csv",
        f"{RESULTS}/differential_expression/PCA_conditions.pdf",
        f"{RESULTS}/differential_expression/volcano.pdf",
        f"{RESULTS}/differential_expression/heatmap.pdf"


# -------------------------
# Fast QC (Pre-Trimming)
# -------------------------
rule fastqc_raw:
    input:
        r1=f"{RAW}/{{sample}}_1.fastq.gz",
	r2=f"{RAW}/{{sample}}_2.fastq.gz"
    output:
        html1=f"{RESULTS}/qc/raw/{{sample}}_1_fastqc.html",
        zip1=f"{RESULTS}/qc/raw/{{sample}}_1_fastqc.zip",
        html2=f"{RESULTS}/qc/raw/{{sample}}_2_fastqc.html",
        zip2=f"{RESULTS}/qc/raw/{{sample}}_2_fastqc.zip"
    threads: 4
    shell:
        """
        mkdir -p {RESULTS}/qc/raw
        fastqc {input.r1} {input.r2} \
        -o {RESULTS}/qc/raw \
        --threads {threads}
        """


# -------------------------
# Trimming 
# -------------------------
rule trim:
    input:
        r1=f"{RAW}/{{sample}}_1.fastq.gz",
        r2=f"{RAW}/{{sample}}_2.fastq.gz"
    output:
        r1=f"{RESULTS}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"{RESULTS}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
        unp1=f"{RESULTS}/trimmed/{{sample}}_R1.unpaired.fastq.gz",
        unp2=f"{RESULTS}/trimmed/{{sample}}_R2.unpaired.fastq.gz"
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
        r1=f"{RESULTS}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"{RESULTS}/trimmed/{{sample}}_R2.trimmed.fastq.gz"
    output:
        html1=f"{RESULTS}/qc/trimmed/{{sample}}_R1.trimmed_fastqc.html",
        html2=f"{RESULTS}/qc/trimmed/{{sample}}_R2.trimmed_fastqc.html",
        zip1=f"{RESULTS}/qc/trimmed/{{sample}}_R1.trimmed_fastqc.zip",
        zip2=f"{RESULTS}/qc/trimmed/{{sample}}_R2.trimmed_fastqc.zip"
    threads: 4
    shell:
        """
        mkdir -p {RESULTS}/qc/trimmed

        fastqc {input.r1} {input.r2} \
            --outdir {RESULTS}/qc/trimmed \
            --threads {threads}
        """

# -------------------------
# HISAT2 index
# -------------------------
rule hisat2_index:
    input:
        config["reference"]["genome"]
    output:
        directory(HISAT2_INDEX)
    threads: 4
    log:
        f"{LOGS}/hisat2_index.log"
    shell:
        """
        mkdir -p {HISAT2_INDEX} {LOGS}

        hisat2-build {input} {HISAT2_INDEX}/mouse_index \
        > {log} 2>&1
        """


# -------------------------
# HISAT2 alignment
# -------------------------
rule align:
    input:
        r1=f"{RESULTS}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"{RESULTS}/trimmed/{{sample}}_R2.trimmed.fastq.gz", 
	index=rules.hisat2_index.output
    output:
        bam=f"{RESULTS}/alignments/{{sample}}.sorted.bam",
	bai=f"{RESULTS}/alignments/{{sample}}.sorted.bam.bai"
    threads: 4
    shell:
        """
        hisat2 \
            -x {HISAT2_INDEX}/mouse_index \
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
        f"{RESULTS}/alignments/{{sample}}.sorted.bam"
    output:
        f"{RESULTS}/qc/alignment/{{sample}}.flagstat.txt"
    log:
        f"{LOGS}/flagstat/{{sample}}.log"
    shell:
        """
        mkdir -p {RESULTS}/qc/alignment {LOGS}/flagstat

        samtools flagstat {input} > {output}

        echo "Done" > {log}
        """

# -------------------------
# featureCounts
# -------------------------
rule featurecounts:
    input:
        bams=expand(
            f"{RESULTS}/alignments/{{sample}}.sorted.bam",
            sample=SAMPLES
        ),
        gtf=config["reference"]["gtf"]
    output:
        f"{RESULTS}/counts/gene_counts.txt",
        f"{RESULTS}/counts/gene_counts.txt.summary"
    threads: 4
    shell:
        """
        mkdir -p {RESULTS}/counts {LOGS}

        featureCounts -T {threads} \
        -p \
        -a {input.gtf} \
        -o {RESULTS}/counts/gene_counts.txt \
        -t exon \
        -g gene_id \
        {input.bams}
        """


# -------------------------
# DESeq2 Analysis
# -------------------------
rule deseq2:
    input:
        counts=f"{RESULTS}/counts/gene_counts.txt",
        metadata=META
    output:
        de_results=f"{RESULTS}/differential_expression/DE_results.csv",
        comparison=f"{RESULTS}/differential_expression/comparison.csv",
        pca=f"{RESULTS}/differential_expression/PCA_conditions.pdf",
        volcano=f"{RESULTS}/differential_expression/volcano.pdf",
        heatmap=f"{RESULTS}/differential_expression/heatmap.pdf"
    params:
        reference_condition=config["deseq2"]["reference_condition"],
        comparison_condition=config["deseq2"]["comparison_condition"]
    log:
        f"{LOGS}/deseq2.log"
    script:
        f"{SCRIPTS}/deseq2_analysis.R"


# -------------------------
# MultiQC
# -------------------------
rule multiqc:
    input:
        expand(
            f"{RESULTS}/qc/raw/{{sample}}_1_fastqc.zip",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS}/qc/raw/{{sample}}_2_fastqc.zip",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS}/qc/trimmed/{{sample}}_R1.trimmed_fastqc.zip",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS}/qc/trimmed/{{sample}}_R2.trimmed_fastqc.zip",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS}/qc/alignment/{{sample}}.flagstat.txt",
            sample=SAMPLES
        ),
        f"{RESULTS}/counts/gene_counts.txt"
    output:
        f"{RESULTS}/multiqc/multiqc_report.html"
    log:
        f"{LOGS}/multiqc.log"
    shell:
        """
        mkdir -p {RESULTS}/multiqc {LOGS}

        multiqc {RESULTS} \
            --outdir {RESULTS}/multiqc \
            --filename multiqc_report.html \
            --force \
            > {log} 2>&1
        """
