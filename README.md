RNAQuest: Mouse RNA-seq Differential Expression Analysis
Overview

RNAQuest is a Snakemake-based pipeline for bulk RNA-seq analysis of mouse transcriptomic data. It processes paired-end sequencing reads through quality control, trimming, genome alignment, gene-level read counting, and differential expression analysis using DESeq2.

The current workflow is designed for a two-condition differential expression comparison. The reference and comparison conditions are defined in config/config.yaml, while sample identifiers and condition assignments are provided in metadata/samples.csv.

The pipeline generates quality-control reports, alignment and count files, differential expression results, and visualizations to support downstream interpretation.

Experimental Design

The current analysis compares two experimental conditions:

Setting	Description
Reference condition	Defined by deseq2.reference_condition in the configuration
Comparison condition	Defined by deseq2.comparison_condition in the configuration
Sample metadata	metadata/samples.csv
Statistical analysis	DESeq2

The number of samples and biological replicates is determined by the metadata file. Ensure that the sample identifiers in the metadata match the corresponding raw FASTQ filenames.

Project Structure
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
    └── ...


This tree represents the intended workflow organization. Additional intermediate files, reference indexes, and test artifacts may be present in the working directory.

Analysis Workflow
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


The exact columns and filtering criteria in the result tables are determined by scripts/deseq2_analysis.R.

7. MultiQC Reporting

MultiQC aggregates quality-control results from raw and trimmed reads, alignment statistics, and gene-counting outputs where supported by the input files.

The main report is:

results/multiqc/multiqc_report.html

Configuration and Metadata
Configuration

The main configuration file is:

config/config.yaml


It defines the input and output paths, reference files, scripts, and DESeq2 comparison settings.

The two-condition comparison is controlled by:

deseq2:
  reference_condition: ...
  comparison_condition: ...


Use the actual condition names defined in the sample metadata. The configuration and metadata must agree on the condition labels.

Sample Metadata

Sample information is maintained in:

metadata/samples.csv


The Snakefile reads sample identifiers from the sample column. Each identifier must correspond to the paired-end raw FASTQ files expected by the workflow:

data/raw/<sample>_1.fastq.gz
data/raw/<sample>_2.fastq.gz


The metadata must also contain the condition information required by the DESeq2 analysis script.

Running the Analysis
1. Set Up the Environment

Create the Conda environment from the provided environment definition:

conda env create -f environment.yml


Activate the environment using its configured name:

conda activate RNAseq


Replace RNAseq if the environment name in environment.yml differs.

2. Run the Complete Workflow

From the project root directory, run:

snakemake --cores 4


Replace 4 with the number of CPU cores available for the analysis.

Snakemake determines which outputs need to be generated based on the workflow dependencies and existing files.

3. Run the CI Test Suite

The repository includes a test suite under tests/, with scripts for creating test data, preparing a miniature reference, running tests, and validating outputs.

Run the test runner from the project root:

bash tests/scripts/run_tests.sh


The test suite is intended to verify the workflow using test data rather than requiring the complete production dataset. Review the test script for its exact prerequisites and validation criteria.

Differential Expression Interpretation

Differential expression results typically include the following DESeq2 statistics, depending on the output columns retained by the analysis script:

baseMean: average normalized expression across samples.
log2FoldChange: estimated log2 fold change between the comparison and reference conditions.
pvalue: unadjusted statistical significance.
padj: multiple-testing-adjusted p-value.

A commonly used significance threshold is:

padj < 0.05
|log2FoldChange| >= 1


These thresholds should only be interpreted as the pipeline's actual filtering criteria if they are implemented in the R analysis script. Otherwise, they can be applied during downstream result interpretation.

Biological replicate numbers, sample quality, and experimental design affect statistical power. Differential expression findings should be interpreted in the context of the study and validated where appropriate.

Continuous Integration and Testing

The tests/ directory provides supporting resources for testing the pipeline, including:

Test configuration and metadata.
Test input and reference data.
Miniature reference preparation scripts.
Synthetic read generation utilities.
Test execution and output validation scripts.

The CI workflow should verify that the two-condition pipeline completes successfully and generates the required outputs, including the count matrix, MultiQC report, differential expression tables, PCA plot, volcano plot, and heatmap.

Reproducibility

Important files for reproducing the analysis include:

Snakefile
config/config.yaml
environment.yml
metadata/samples.csv
scripts/deseq2_analysis.R
tests/

To record the Conda environment:

conda env export > environment_full.yml


To record R session information:

mkdir -p logs
Rscript -e "sessionInfo()" > logs/R_sessionInfo.txt

Notes
The workflow is configured for mouse bulk RNA-seq data and paired-end sequencing reads.
Reference genome, annotation, and sample paths must match the configuration.
Intermediate files and logs are useful for troubleshooting and reproducibility.
Test data and production data should be kept distinct.
Final biological interpretation may be supplemented with pathway enrichment analysis and independent validation of candidate genes.
