# RNAQuest: Mouse RNA-seq Differential Expression Analysis

## Overview

RNAQuest is a Snakemake-based pipeline for bulk RNA-seq analysis of mouse transcriptomic data. It processes paired-end sequencing reads through quality control, trimming, genome alignment, gene-level read counting, and differential expression analysis using DESeq2.

The current workflow is designed for a **two-condition differential expression comparison**. The reference and comparison conditions are defined in `config/config.yaml`, while sample identifiers and condition assignments are provided in `metadata/samples.csv`.

The pipeline generates quality-control reports, alignment and count files, differential expression results, and visualizations to support downstream interpretation.

## Experimental Design

The current analysis compares two experimental conditions:

| Setting | Description |
|---|---|
| Reference condition | Defined by `deseq2.reference_condition` in the configuration |
| Comparison condition | Defined by `deseq2.comparison_condition` in the configuration |
| Sample metadata | `metadata/samples.csv` |
| Statistical analysis | DESeq2 |

The number of samples and biological replicates is determined by the metadata file. Ensure that the sample identifiers in the metadata match the corresponding raw FASTQ filenames.

## Project Structure

```text
RNAQuest/
├── Snakefile
├── README.md
├── config/
│   └── config.yaml
├── environment.yml
├── metadata/
│   └── samples.csv
├── scripts/
│   └── deseq2_analysis.R
├── data/
│   ├── raw/
│   ├── trimmed/
│   ├── reference/
│   └── adapters/
├── results/
│   ├── qc/
│   │   ├── raw/
│   │   ├── trimmed/
│   │   └── alignment/
│   ├── trimmed/
│   ├── alignments/
│   ├── counts/
│   │   └── gene_counts.txt
│   ├── differential_expression/
│   │   ├── DE_results.csv
│   │   ├── comparison.csv
│   │   ├── PCA_conditions.pdf
│   │   ├── volcano.pdf
│   │   └── heatmap.pdf
│   └── multiqc/
│       └── multiqc_report.html
├── logs/
└── tests/
    ├── config/
    ├── data/
    ├── metadata/
    ├── results/
    ├── scripts/
    │   ├── create_mini_reference.sh
    │   ├── create_test_dataset.sh
    │   ├── generate_synthetic_reads.R
    │   ├── run_tests.sh
    │   └── validate.py
    └── workflows
```

This tree represents the intended workflow organization. Additional intermediate files, reference indexes, and test artifacts may be present in the working directory.

## Analysis Workflow
1. Raw Read Quality Control

FastQC evaluates the raw paired-end FASTQ files before trimming.

Outputs are written to:

results/qc/raw/

2. Read Trimming

Trimmomatic removes low-quality regions and trims reads using the following settings:

Sliding-window trimming: SLIDINGWINDOW:4:20
Minimum read length: MINLEN:50
Paired-end processing with unpaired reads retained

Outputs are written to:

results/trimmed/

3. Post-Trimming Quality Control

FastQC evaluates the trimmed paired-end reads to assess read quality after preprocessing.

Outputs are written to:

results/qc/trimmed/

4. Genome Alignment

HISAT2 aligns trimmed reads to the configured mouse reference genome. SAMtools sorts and indexes the resulting BAM files.

Reference genome and annotation paths are defined in config/config.yaml.

The genome reference uses the mouse mm39 assembly when configured accordingly.

Outputs include:

results/alignments/


Sorted BAM files and their corresponding BAM indexes are generated for each sample. Alignment statistics are collected using samtools flagstat and stored in:

results/qc/alignment/

5. Gene-Level Read Counting

featureCounts summarizes reads overlapping annotated exons and groups counts by gene_id.

The annotation file is supplied through the pipeline configuration.

Main outputs:

results/counts/gene_counts.txt
results/counts/gene_counts.txt.summary

6. Differential Expression Analysis

DESeq2 performs differential expression analysis using the two conditions specified in config/config.yaml.

The analysis uses the gene count matrix and sample metadata to estimate expression differences between the reference and comparison conditions.

Outputs:

File	Description
DE_results.csv	Differential expression results
comparison.csv	Comparison summary or associated analysis output
PCA_conditions.pdf	Principal component analysis visualization
volcano.pdf	Volcano plot of differential expression
heatmap.pdf	Heatmap visualization

All files are written to:

results/differential_expression/
