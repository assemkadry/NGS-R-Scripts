# DESeq2 Differential Expression Pipeline

This repository contains a minimal, streamlined pipeline for differential expression analysis using RNA-seq count data. It includes:

- Differential expression analysis with **DESeq2**
- Visualization with **Volcano Plot**
- Visualization with **Heatmap**

---

## 📂 Files

| Script            | Description |
|------------------|-------------|
| `DESeq2.r`        | Runs DESeq2 on a count matrix (via stdin). Outputs a results table and normalized counts for significant DEGs. |
| `volcano_plot.r`  | Creates a volcano plot PDF from DESeq2 results (via stdin). Highlights up/down-regulated genes. |
| `heatmap.r`       | Generates a z-score normalized heatmap (via stdin) of differentially expressed genes. Outputs a PDF. |

---

## 🔧 Requirements

These R packages must be installed:

```r
install.packages(c("ggplot2", "gplots", "RColorBrewer"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
```

---

## ▶️ Usage

### 1. Run DESeq2

```bash
cat counts.tsv | Rscript DESeq2.r 3x3
```

- Assumes 3 replicates per group (NxM format).
- Outputs:
  - `results_deseq2.tsv`
  - `norm-matrix-deseq2.tsv`

---

### 2. Generate a Volcano Plot

```bash
cat results_deseq2.tsv | Rscript volcano_plot.r > volcano_plot.pdf
```

- Uses adjusted p-values (`padj`) and log2 fold change.
- Highlights significant genes.

---

### 3. Generate a Heatmap

```bash
cat norm-matrix-deseq2.tsv | Rscript heatmap.r > heatmap.pdf
```

- Z-score scaled across rows (genes).
- Samples are clustered and visualized.

---

## 📝 Notes

- All scripts accept input via `stdin` for easy shell integration.
- Sample column names should follow a recognizable pattern (e.g., `cancer_1`, `normal_1`).
- Output is always written to standard output or named files.

---

## 🧪 Sample Input Format

**counts.tsv**
```
gene	sample1	sample2	sample3	sample4	sample5	sample6
GENE1	10	20	12	200	210	190
GENE2	0	1	0	150	140	160
...
```

---

## 📄 License

This project is licensed under the MIT License.
