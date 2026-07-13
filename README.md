# RNAQuest: Mouse RNA-seq Differential Expression Analysis

## Overview

This project contains an RNA-seq analysis pipeline for mouse transcriptomic data. The workflow processes raw sequencing reads through quality control, trimming, alignment, read counting, and differential expression analysis using DESeq2.

The primary goal of this analysis is to identify genes with altered expression between experimental conditions:

* Sed: Sedentary control
* Ex: Exercise
* ExSed: Exercise + Sedentary condition

## Project Structure

```
RNAQuest/
├── Snakefile                  # Snakemake workflow
├── config/                    # Pipeline configuration files
├── environment.yml            # Software environment
├── scripts/                   # Analysis scripts
│   └── deseq2_analysis.R      # DESeq2 differential expression script
├── metadata/
│   └── samples.csv            # Sample information and conditions
├── data/
│   ├── raw/                   # Raw sequencing reads
│   ├── trimmed/               # Trimmed reads
│   ├── reference/             # Genome reference and annotation
│   └── adapters/              # Adapter sequences
├── results/
│   ├── alignments/            # Aligned BAM/SAM files
│   ├── counts/                # featureCounts output
│   ├── differential_expression/
│   │   ├── Ex_vs_Sed.csv
│   │   ├── ExSed_vs_Sed.csv
│   │   ├── ExSed_vs_Ex.csv
│   │   ├── PCA_conditions.pdf
│   │   └── volcano_Ex_vs_Sed.pdf
│   └── multiqc/               # MultiQC reports
└── logs/                      # Pipeline execution logs
```

## Experimental Design

Samples were analyzed in three groups:

| Condition | Replicates |
| --------- | ---------- |
| Sed       | 2          |
| Ex        | 2          |
| ExSed     | 2          |

The differential expression analysis used DESeq2 with condition as the experimental factor.

## Analysis Workflow

### 1. Quality Control

Raw sequencing reads were evaluated using FastQC.

### 2. Read Trimming

Adapters and low-quality regions were removed before alignment.

### 3. Alignment

Reads were aligned to the mouse reference genome using the configured aligner.

Reference files:

* Genome: mouse mm39
* Annotation: GTF file in `data/reference/`

### 4. Read Counting

Gene-level counts were generated using featureCounts.

Output:

```
results/counts/gene_counts.txt
```

### 5. Differential Expression Analysis

DESeq2 was used to compare:

* Exercise vs Sedentary

```
results/differential_expression/Ex_vs_Sed.csv
```

* Exercise + Sedentary vs Sedentary

```
results/differential_expression/ExSed_vs_Sed.csv
```

* Exercise + Sedentary vs Exercise

```
results/differential_expression/ExSed_vs_Ex.csv
```

## Running the Analysis

Activate the conda environment:

```bash
conda activate RNAseq
```

Run the complete workflow:

```bash
snakemake --cores <number_of_cores>
```

Run the DESeq2 analysis manually:

```bash
Rscript scripts/deseq2_analysis.R
```

## Differential Expression Interpretation

Results contain:

* `baseMean`: average normalized expression
* `log2FoldChange`: expression change between conditions
* `pvalue`: statistical significance before correction
* `padj`: multiple-testing corrected significance value

Significance criteria commonly used:

```
padj < 0.05
|log2FoldChange| >= 1
```

Because this dataset contains two biological replicates per condition, results should be interpreted carefully. Genes with strong fold changes but higher adjusted p-values may represent candidate genes requiring further validation.

## Results Summary

Generated outputs include:

* Gene count matrix
* Differential expression tables
* PCA visualization
* Volcano plots
* MultiQC quality reports

## Reproducibility

Important files for reproducing this analysis:

* `Snakefile`
* `config/`
* `environment.yml`
* `metadata/samples.csv`
* `scripts/deseq2_analysis.R`

Software versions and package information can be recorded using:

```bash
conda env export > environment_full.yml
```

and:

```bash
Rscript -e "sessionInfo()" > logs/R_sessionInfo.txt
```

## Notes

* All analysis was performed on Ubuntu Linux.
* Intermediate files are retained to allow troubleshooting and pipeline verification.
* Final biological interpretation should include pathway enrichment analysis and validation of candidate genes.
