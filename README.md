# RNAQuest

RNAQuest is a bulk RNA-seq analysis pipeline for mouse transcriptomic data, built with Snakemake. It processes paired-end FASTQ files through quality control, read trimming, transcript-level quantification with Salmon, and gene-level differential expression analysis using DESeq2 and tximport.

The workflow is configurable through YAML configuration files and sample metadata, with automated quality-control reporting and differential expression visualizations.

The pipeline is designed to compare two experimental conditions, with the reference and comparison conditions specified in the configuration.
# Features

- Automated RNA-seq analysis using Snakemake
- Raw read quality control with FastQC
- Paired-end read trimming with Trimmomatic
- Post-trimming quality control with FastQC
- Transcript-level quantification using Salmon
- Lightweight transcript quantification with selective alignment-based mapping validation
- Transcript-to-gene summarization using tximport
- Two-condition differential expression analysis using DESeq2
- Principal component analysis (PCA), volcano plot, and heatmap generation
- Aggregated quality-control reporting with MultiQC
- Configurable sample metadata, transcript references, Salmon indexes, and transcript-to-gene mappings
- CI testing support using synthetic data and a miniature reference
- Reproducible Conda-based software environment

# Tools Used

| Step | Tool |
|---|---|
| Workflow Management | Snakemake |
| Quality Control | FastQC |
| Read Trimming | Trimmomatic |
| Transcript Quantification | Salmon |
| Text-to-Gene Summarization | tximport |
| Differential Expression | DESeq2 |
| QC Aggregation | MultiQC |
| Statistical Analysis and Visualization | R |

# Pipeline Workflow

```text
Paired-end FASTQ files
          |
          v
   FastQC (raw reads)
          |
          v
     Trimmomatic
          |
          v
 FastQC (trimmed reads)
          |
          v
   Salmon quantification
          |
          v
 Transcript-level quant.sf
          |
          v
 tximport: transcript-to-gene
       summarization
          |
          v
 DESeq2 (two-condition comparison)
          |
          +-------------------+
          |         |         |
          v         v         v
         PCA     Volcano    Heatmap

FastQC + Salmon outputs
          |
          v
       MultiQC report

```


## Project Structure

```text
RNAQuest/
├── config/
│   └── config.yaml
├── data/
│   ├── adapters/
│   ├── raw/
│   ├── reference/
│   └── trimmed/
├── metadata/
│   └── samples.csv
├── results/
│   ├── qc/
│   │   ├── raw/
│   │   └── trimmed/
│   ├── trimmed/
│   ├── salmon/
│   │   └── <sample>/
│   │       ├── quant.sf
│   │       └── aux_info/
│   │           └── meta_info.json
│   ├── differential_expression/
│   │   ├── DE_results.csv
│   │   ├── comparison.csv
│   │   ├── PCA_conditions.pdf
│   │   ├── volcano.pdf
│   │   └── heatmap.pdf
│   └── multiqc/
│       └── multiqc_report.html
├── logs/
│   ├── salmon_index.log
│   ├── salmon_quant/
│   ├── deseq2.log
│   └── multiqc.log
├── scripts/
│   └── deseq2_analysis.R
├── tests/
│   ├── config/
│   │   └── config.yaml
│   ├── data/
│   │   ├── metadata.csv
│   │   ├── raw/
│   │   └── reference/
│   ├── metadata/
│   │   └── samples.csv
│   ├── results/
│   └── scripts/
│       ├── create_mini_reference.sh
│       ├── create_test_dataset.sh
│       ├── filter_test_gtf.py
│       ├── filter_transcripts.py
│       ├── generate_synthetic_reads.R
│       ├── deseq2_analysis.R
│       ├── run_tests.sh
│       └── validate.py
├── Snakefile
├── environment.yml
├── LICENSE
└── README.md
```
Unlike the previous HISAT2-based workflow, the current pipeline uses Salmon to quantify transcript abundance without producing genome-aligned BAM files. Transcript-level estimates are summarized to the gene level using tximport before differential expression analysis with DESeq2.

## Installation
Clone the repository and enter the project directory:

```bash
git clone https://github.com/yourusername/RNAQuest.git
cd RNAQuest
```
Create the Conda environment:

```bash
conda env create -f environment.yml
```
Activate the environment:

```bash
conda activate RNAseq
```
Replace RNAseq if the environment name defined in environment.yml differs.

## Configuration
The main configuration file is:

```text
config/config.yaml
```
It defines the input and output paths, reference genome, annotation file, script paths, and DESeq2 comparison settings.

The two-condition comparison is configured using:

```yaml
deseq2:
  reference_condition: ...
  comparison_condition: ...
```
Replace the placeholders with the actual condition labels used in the sample metadata.

Sample information is maintained in:

```text
metadata/samples.csv
```
The Snakefile reads sample identifiers from the sample column. Each sample must have matching paired-end FASTQ files following this naming convention:

```text
data/raw/<sample>_1.fastq.gz
data/raw/<sample>_2.fastq.gz
```
Ensure that the sample identifiers and condition labels agree with the configuration and the requirements of the DESeq2 analysis script.

Run Pipeline
Production Analysis
Run the complete workflow from the project root:

```bash
snakemake --cores 4
```
Replace 4 with the number of CPU cores available.

Snakemake uses the configured sample metadata and workflow dependencies to determine which outputs need to be generated.

## License

MIT License

