#!/usr/bin/env Rscript

suppressMessages({
  library(gplots)
  library(RColorBrewer)
})

# Read from stdin
data <- read.table("stdin", header = TRUE, sep = "\t", row.names = 1)

# Matrix of expression values
mat <- as.matrix(data)

# Optional: add small noise to avoid 0 variance
mat <- jitter(mat, factor = 1, amount = 1e-5)

# Z-score normalization by gene (row)
zscore <- t(scale(t(mat)))

# Set colors
colors <- colorRampPalette(c("blue", "white", "red"))(256)

# Plot to PDF via stdout pipe
pdf("|cat", width = 8, height = 10)
heatmap.2(zscore,
          col = colors,
          trace = "none",
          density.info = "none",
          margins = c(8, 8),
          cexRow = 0.6,
          cexCol = 0.8,
          srtCol = 45,
          scale = "none",
          key = TRUE,
          key.title = "Z-score",
          key.xlab = "Expression")
dev.off()