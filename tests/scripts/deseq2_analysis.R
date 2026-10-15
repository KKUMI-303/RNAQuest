suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
})

# Inputs from Snakemake

tx2gene_file <- snakemake@input[["tx2gene"]]
metadata_file <- snakemake@input[["metadata"]]

reference_condition <- snakemake@params[["reference_condition"]]
comparison_condition <- snakemake@params[["comparison_condition"]]

de_results_file <- snakemake@output[["de_results"]]
comparison_file <- snakemake@output[["comparison"]]
pca_file <- snakemake@output[["pca"]]
volcano_file <- snakemake@output[["volcano"]]
heatmap_file <- snakemake@output[["heatmap"]]

quant_files <- unlist(snakemake@input[["quant"]])

sample_ids <- strsplit(
  snakemake@params[["samples"]],
  ",",
  fixed = TRUE
)[[1]]

if (length(quant_files) != length(sample_ids)) {
  stop("Number of Salmon quantification files does not match sample IDs.")
}

names(quant_files) <- sample_ids

de_dir <- dirname(de_results_file)
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

# Validate Salmon quantification files
if (length(quant_files) == 0 || any(!file.exists(quant_files))) {
  stop("One or more Salmon quant.sf files are missing.")
}

# Load metadata
coldata <- read.csv(
  metadata_file,
  row.names = 1,
  check.names = FALSE
)

if (anyDuplicated(rownames(coldata))) {
  stop("Sample IDs in metadata must be unique.")
}

if (!all(names(quant_files) %in% rownames(coldata))) {
  stop("A Salmon sample ID is missing from the metadata.")
}

coldata <- coldata[names(quant_files), , drop = FALSE]

if (!"condition" %in% colnames(coldata)) {
  stop("Metadata must contain a condition column.")
}

coldata$condition <- factor(
  coldata$condition,
  levels = c(reference_condition, comparison_condition)
)

if (anyNA(coldata$condition)) {
  stop("Metadata contains unexpected or missing condition labels.")
}

if (any(table(coldata$condition) < 3)) {
  stop("Expected at least 3 biological replicates per condition.")
}

# Load transcript-to-gene mapping
tx2gene <- read.csv(
  tx2gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("transcript_id", "gene_id") %in% colnames(tx2gene))) {
  stop("tx2gene.csv must contain transcript_id and gene_id columns.")
}

tx2gene <- unique(tx2gene[, c("transcript_id", "gene_id")])

if (
  anyNA(tx2gene) ||
  any(tx2gene$transcript_id == "") ||
  any(tx2gene$gene_id == "")
) {
  stop("tx2gene.csv contains missing transcript or gene IDs.")
}

# Import Salmon estimates and aggregate to genes
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no",
  dropInfReps = TRUE
)

# Construct DESeq2 dataset
dds <- DESeqDataSetFromTximport(
  txi,
  colData = coldata,
  design = ~ condition
)

# Remove genes with zero counts across all samples
dds <- dds[rowSums(counts(dds)) > 0, ]

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(
    "condition",
    comparison_condition,
    reference_condition
  )
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

write.csv(res_df, de_results_file, row.names = FALSE)
write.csv(res_df, comparison_file, row.names = FALSE)

# Variance-stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# PCA plot
p <- plotPCA(vsd, intgroup = "condition")

ggsave(
  pca_file,
  p,
  width = 7,
  height = 5
)

# Volcano plot
plot_df <- res_df[
  !is.na(res_df$log2FoldChange) & !is.na(res_df$pvalue),
  ,
  drop = FALSE
]

pdf(volcano_file, width = 7, height = 5)

print(
  ggplot(plot_df, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(
      x = "Log2 fold change",
      y = "-Log10 p-value",
      title = paste(comparison_condition, "vs", reference_condition)
    )
)

dev.off()

# Heatmap of the most variable genes
sample_order <- order(coldata$condition)

mat <- assay(vsd)[, sample_order, drop = FALSE]

rv <- matrixStats::rowVars(mat)

selected <- order(rv, decreasing = TRUE)[
  seq_len(min(500, length(rv)))
]

mat <- mat[selected, , drop = FALSE]

# Remove constant rows before scaling
row_sd <- apply(mat, 1, sd)
mat <- mat[row_sd > 0, , drop = FALSE]
mat <- t(scale(t(mat)))

annotation_col <- data.frame(
  Condition = coldata$condition[sample_order]
)

rownames(annotation_col) <- colnames(mat)

pheatmap(
  mat,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  filename = heatmap_file,
  width = 8,
  height = 10
)
