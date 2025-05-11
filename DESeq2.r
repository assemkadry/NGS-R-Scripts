#!/usr/bin/env Rscript

# Usage: cat counts.txt | Rscript deseq2_clean.R 3x3

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript deseq2_clean.R NxM (e.g., 3x3)", call. = FALSE)
}

# Parse design
design <- strsplit(args[1], "x")[[1]]
cond1_num <- as.integer(design[1])
cond2_num <- as.integer(design[2])
conditions <- factor(c(rep("cond1", cond1_num), rep("cond2", cond2_num)))

suppressMessages(library(DESeq2))

# Read count data
count_data <- read.table(file("stdin"), header = TRUE, row.names = 1, sep = "\t")
count_matrix <- round(as.matrix(count_data))
mode(count_matrix) <- "integer"

# Metadata
col_data <- data.frame(condition = conditions)
rownames(col_data) <- colnames(count_matrix)

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

# Results
res <- results(dds, contrast = c("condition", "cond1", "cond2"))
res <- res[order(res$padj), ]

# Format numbers safely
res_df <- as.data.frame(res)
res_df <- format(res_df, digits = 6, scientific = TRUE)
res_df <- cbind(id = rownames(res_df), res_df)

# Save results safely
write.table(res_df, file = "results_deseq2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Save normalized counts of DE genes
de_ids <- rownames(subset(res, padj < 0.05 & !is.na(padj)))
norm_counts <- counts(dds, normalized = TRUE)
norm_df <- data.frame(id = rownames(norm_counts), norm_counts)
de_df <- norm_df[norm_df$id %in% de_ids, ]
write.table(de_df, file = "norm-matrix-deseq2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
