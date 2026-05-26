library(DESeq2)

counts <- read.delim("results/counts/gene_counts.txt", comment.char="#", row.names=1)

counts <- counts[,6:ncol(counts)]  # remove metadata columns

coldata <- data.frame(
  row.names = colnames(counts),
  condition = c("control","treated")
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res),
          "results/differential_expression/DE_results.csv")

# PCA
pdf("results/differential_expression/PCA.pdf")
plotPCA(vst(dds), intgroup="condition")
dev.off()

# Volcano plot
library(ggplot2)
res_df <- as.data.frame(res)

pdf("results/differential_expression/volcano.pdf")
ggplot(res_df, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point() +
  theme_minimal()
dev.off()
