configfile: "config/config.yaml"

import csv
from snakemake.io import directory

RAW = config["paths"]["raw"]
RESULTS = config["paths"]["results"]
LOGS = config["paths"]["logs"]
META = config["paths"]["metadata"]
SCRIPTS = config["scripts"]["test_script"]

TRANSCRIPTS = config["reference"]["transcripts"]
SALMON_INDEX = config["reference"]["salmon_index"]
TX2GENE = config["reference"]["tx2gene"]

with open(META) as f:
    reader = csv.DictReader(f)
    SAMPLES = [row["sample"] for row in reader]

rule all:
    input:
        expand(
            f"{RESULTS}/salmon/{{sample}}/quant.sf",
            sample=SAMPLES
        ),
        f"{RESULTS}/multiqc/multiqc_report.html",
        f"{RESULTS}/differential_expression/DE_results.csv",
        f"{RESULTS}/differential_expression/comparison.csv",
        f"{RESULTS}/differential_expression/PCA_conditions.pdf",
        f"{RESULTS}/differential_expression/volcano.pdf",
        f"{RESULTS}/differential_expression/heatmap.pdf"


# -------------------------
# FastQC: Raw reads
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
    threads: 2
    shell:
        """
        mkdir -p {RESULTS}/qc/raw

        fastqc {input.r1} {input.r2} \
            --outdir {RESULTS}/qc/raw \
            --threads {threads}
        """


# -------------------------
# Paired-end trimming
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
    threads: 2
    shell:
        """
        mkdir -p {RESULTS}/trimmed

        trimmomatic PE \
            -threads {threads} \
            -phred33 \
            {input.r1} {input.r2} \
            {output.r1} {output.unp1} \
            {output.r2} {output.unp2} \
            SLIDINGWINDOW:4:20 MINLEN:50
        """


# -------------------------
# FastQC: Trimmed reads
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
    threads: 2
    shell:
        """
        mkdir -p {RESULTS}/qc/trimmed

        fastqc {input.r1} {input.r2} \
            --outdir {RESULTS}/qc/trimmed \
            --threads {threads}
        """


# -------------------------
# Salmon transcriptome index
# -------------------------
rule salmon_index:
    input:
        transcripts=TRANSCRIPTS
    output:
        index=directory(SALMON_INDEX)
    log:
        f"{LOGS}/salmon_index.log"
    threads: 2
    shell:
        """
        mkdir -p {LOGS}

        salmon index \
            -t {input.transcripts} \
            -i {output.index} \
            -k 21 \
            > {log} 2>&1
        """


# -------------------------
# Salmon paired-end quantification
# -------------------------
rule salmon_quant:
    input:
        r1=f"{RESULTS}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2=f"{RESULTS}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
        index=rules.salmon_index.output
    output:
        quant=f"{RESULTS}/salmon/{{sample}}/quant.sf",
        meta=f"{RESULTS}/salmon/{{sample}}/aux_info/meta_info.json"
    log:
        f"{LOGS}/salmon_quant/{{sample}}.log"
    threads: 2
    shell:
        """
        mkdir -p {RESULTS}/salmon/{wildcards.sample}
        mkdir -p {LOGS}/salmon_quant

        salmon quant \
            -i {input.index} \
            -l A \
            -1 {input.r1} \
            -2 {input.r2} \
            --validateMappings \
            -p {threads} \
            -o {RESULTS}/salmon/{wildcards.sample} \
            > {log} 2>&1
        """


# -------------------------
# DESeq2 using tximport
# -------------------------

rule deseq2:
    input:
        quant=expand(
            f"{RESULTS}/salmon/{{sample}}/quant.sf",
            sample=SAMPLES
        ),
        tx2gene=TX2GENE,
        metadata=META
    output:
        de_results=f"{RESULTS}/differential_expression/DE_results.csv",
        comparison=f"{RESULTS}/differential_expression/comparison.csv",
        pca=f"{RESULTS}/differential_expression/PCA_conditions.pdf",
        volcano=f"{RESULTS}/differential_expression/volcano.pdf",
        heatmap=f"{RESULTS}/differential_expression/heatmap.pdf"
    params:
        samples=",".join(SAMPLES),
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
            f"{RESULTS}/salmon/{{sample}}/quant.sf",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS}/salmon/{{sample}}/aux_info/meta_info.json",
            sample=SAMPLES
        )
    output:
        f"{RESULTS}/multiqc/multiqc_report.html"
    log:
        f"{LOGS}/multiqc.log"
    shell:
        """
        mkdir -p {RESULTS}/multiqc {LOGS}

        multiqc \
            {RESULTS}/qc/raw \
            {RESULTS}/qc/trimmed \
            {RESULTS}/salmon \
            --outdir {RESULTS}/multiqc \
            --filename multiqc_report.html \
            --force \
            > {log} 2>&1
        """
