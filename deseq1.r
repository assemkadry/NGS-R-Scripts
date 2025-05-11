#!/usr/bin/env Rscript

# Run with: cat counts.txt | Rscript deseq2_pipe.R 3x3

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript deseq2_pipe.R NxM  (e.g., 3x3)", call. = FALSE)
}

# Parse design
design <- strsplit(args[1], "x")[[1]]
cond1_num <- as.integer(design[1])
cond2_num <- as.integer(design[2])
conditions <- factor(c(rep("cond1", cond1_num), rep("cond2", cond2_num)))

# Load libraries
suppressMessages({
  library(DESeq2)
})

# Read from stdin
count_data <- read.table(file("stdin"), header = TRUE, row.names = 1, sep = "\t")

# Ensure integer counts
count_matrix <- round(as.matrix(count_data))
mode(count_matrix) <- "integer"

# Create metadata
col_data <- data.frame(condition = conditions)
rownames(col_data) <- colnames(count_matrix)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) > 10, ]  # filter low counts
dds <- DESeq(dds)

# Differential expression
res <- results(dds, contrast = c("condition", "cond1", "cond2"))
res <- res[order(res$padj), ]
write.table(res, file = "", sep = "\t", quote = FALSE)

# Get DEGs with padj < 0.05
de_ids <- rownames(subset(res, padj < 0.05))

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
norm_df <- data.frame(id = rownames(norm_counts), norm_counts)
de_df <- norm_df[norm_df$id %in% de_ids, ]

# Save
write.table(de_df, file = "norm-matrix-deseq2.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
