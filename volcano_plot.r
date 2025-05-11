#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))

# Read from stdin
res <- read.table("stdin", header = TRUE, sep = "\t", row.names = 1)

# Convert columns to numeric
res$log2FoldChange <- as.numeric(res$log2FoldChange)
res$padj <- as.numeric(res$padj)

# Filter out rows with NA padj
res <- res[!is.na(res$padj), ]

# Prepare data
res$gene <- rownames(res)
res$threshold <- "Not Sig"
res$threshold[res$padj < 0.05 & res$log2FoldChange > 1] <- "Up"
res$threshold[res$padj < 0.05 & res$log2FoldChange < -1] <- "Down"

# Cap extreme -log10(padj) for visual clarity
res$log10padj <- -log10(res$padj)
res$log10padj[res$log10padj > 300] <- 300

# Volcano plot
p <- ggplot(res, aes(x = log2FoldChange, y = log10padj, color = threshold)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c("Up" = "firebrick", "Down" = "royalblue", "Not Sig" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", color = "Group")

# Output to PDF via pipe
pdf("|cat", width = 7, height = 5)
print(p)
dev.off()
