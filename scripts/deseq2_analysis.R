library(DESeq2)
library(ggplot2)
library(pheatmap)
library(matrixStats)


counts_file <- snakemake@input[["counts"]]
metadata_file <- snakemake@input[["metadata"]]
reference_condition <- snakemake@params[["reference_condition"]]
comparison_condition <- snakemake@params[["comparison_condition"]]

de_dir <- dirname(snakemake@output[["de_results"]]) 

dir.create(
  de_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Load counts
counts <- read.delim(
  counts_file,
  comment.char = "#",
  row.names = 1,
  check.names = FALSE
)

counts <- counts[, !(colnames(counts) %in%
  c("Chr", "Start", "End", "Strand", "Length"))]

colnames(counts) <- gsub("^.*/", "", colnames(counts))
colnames(counts) <- sub("\\.sorted\\.bam$", "", colnames(counts))
counts <- as.matrix(counts)

# Load metadata
coldata <- read.csv(
  metadata_file,
  row.names = 1,
  check.names = FALSE
)

coldata <- coldata[colnames(counts), , drop = FALSE]

stopifnot(identical(colnames(counts), rownames(coldata)))

# Two-condition design
coldata$condition <- factor(
  coldata$condition,
  levels = c(reference_condition, comparison_condition)
)

if (any(is.na(coldata$condition))) {
  stop("Expected conditions: Control and Treatment")
}

if (any(table(coldata$condition) < 3)) {
  stop("Expected 3 biological replicates per condition")
}

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", comparison_condition, reference_condition)
)

write.csv(
  as.data.frame(res),
  file.path(de_dir, "DE_results.csv")
)

write.csv(
  as.data.frame(res),
  file.path(de_dir, "comparison.csv")
)

# VST
vsd <- vst(dds, blind = FALSE)

# PCA
p <- plotPCA(vsd, intgroup = "condition")

ggsave(
  file.path(de_dir, "PCA_conditions.pdf"),
  p,
  width = 7,
  height = 5
)

# Volcano
res_df <- na.omit(as.data.frame(res))

pdf(file.path(de_dir, "volcano.pdf"))

print(
  ggplot(
    res_df,
    aes(log2FoldChange, -log10(pvalue))
  ) +
    geom_point(alpha = 0.5) +
    theme_minimal()
)

dev.off()

# Heatmap
sample_order <- order(coldata$condition)

mat <- assay(vsd)[, sample_order, drop = FALSE]

rv <- rowVars(mat)
select <- order(rv, decreasing = TRUE)[
  seq_len(min(500, length(rv)))
]

mat <- mat[select, , drop = FALSE]
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
  filename = file.path(de_dir, "heatmap.pdf"),
  width = 8,
  height = 10
)
