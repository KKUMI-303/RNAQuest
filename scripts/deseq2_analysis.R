library(DESeq2)
library(ggplot2)

# -----------------------------
# 1. Load featureCounts output
# -----------------------------
counts <- read.delim(
  "results/counts/gene_counts.txt",
  comment.char = "#",
  row.names = 1,
  check.names = FALSE
)

# REMOVE annotation columns (critical step you removed accidentally)
counts <- counts[, !(colnames(counts) %in%
  c("Chr", "Start", "End", "Strand", "Length"))]

# Remove annotation columns (featureCounts standard format)
colnames(counts) <- gsub("^.*/", "", colnames(counts)) 
colnames(counts) <- sub("\\.sorted\\.bam$", "", colnames(counts)) 
counts <- as.matrix(counts)

# -----------------------------
# 2. Load sample metadata
# -----------------------------
coldata <- read.csv(
  "metadata/samples.csv",
  row.names = 1
)

# enforce identical ordering based on counts
coldata <- coldata[colnames(counts), , drop=FALSE]

# hard validation
stopifnot(setequal(colnames(counts), rownames(coldata)))
stopifnot(all(colnames(counts) == rownames(coldata)))

coldata$condition <- factor(coldata$condition)



# -----------------------------
# 3. Check replicates (critical)
# -----------------------------
tab <- table(coldata$condition)

cat("\nSample distribution:\n")
print(tab)
print(table(coldata$condition))
print(colnames(counts))

if (any(tab < 2)) {
  cat("\n⚠ WARNING: DESeq2 requires at least 2 replicates per condition.\n")
  cat("Switching to simple log2 fold-change analysis only.\n")

  # -----------------------------
  # FALLBACK MODE: no replicates
  # -----------------------------
  res <- data.frame(
    gene = rownames(counts),
    log2FC = log2((counts[,2] + 1) / (counts[,1] + 1))
  )

  write.csv(res,
            "results/differential_expression/DE_results.csv",
            row.names = FALSE)

  # simple PCA (exploratory only)
  pdf("results/differential_expression/PCA.pdf")
  pca <- prcomp(t(log2(counts + 1)))
  plot(pca$x[,1], pca$x[,2],
       col = as.factor(coldata$condition),
       pch = 19,
       xlab = "PC1", ylab = "PC2")
  dev.off()

} else {

  # -----------------------------
  # DESEQ2 MODE (valid case)
  # -----------------------------
  # reorder first (critical)
  coldata <- coldata[colnames(counts), , drop = FALSE]

  # set reference
  coldata$condition <- relevel(
  factor(coldata$condition),
  ref = "Sed"
  )

  dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
  )

  dds <- DESeq(dds)

  # explicit contrasts
  res_ex_vs_sed <- results(dds, contrast=c("condition","Ex","Sed"))
  res_exsed_vs_sed <- results(dds, contrast=c("condition","ExSed","Sed"))
  res_exsed_vs_ex <- results(dds, contrast=c("condition","ExSed","Ex"))

  # save
  write.csv(as.data.frame(res_ex_vs_sed),
          "results/differential_expression/Ex_vs_Sed.csv")

  write.csv(as.data.frame(res_exsed_vs_sed),
          "results/differential_expression/ExSed_vs_Sed.csv")

  write.csv(as.data.frame(res_exsed_vs_ex),
          "results/differential_expression/ExSed_vs_Ex.csv")

  # PCA (valid only here)
  vsd <- vst(dds)

  pdf("results/differential_expression/PCA_conditions.pdf")
  plotPCA(vsd, intgroup = "condition")
  dev.off()

  # Volcano plot (valid only here)
  res_df <- as.data.frame(res_ex_vs_sed)
  res_df <- na.omit(res_df)

  pdf("results/differential_expression/volcano_Ex_vs_Sed.pdf")
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(alpha = 0.5) +
    theme_minimal()
  dev.off()
}