
# RNA-Seq Differential Expression Analysis Pipeline

This repository contains scripts for performing differential gene expression analysis using RNA-seq data. It includes steps for data preprocessing, DESeq2 analysis, visualization (heatmaps and volcano plots), and a shell script for basic environment setup.

## Folder Structure

```
.
├── DESeq2.r           # Differential expression analysis using DESeq2
├── heatmap.r          # Script to generate clustered heatmaps
├── volcano_plot.r     # Script to generate volcano plots
├── RNA_script.sh      # Shell script for data preparation or environment setup
└── README.md          # This file
```

## Requirements

- **R (version ≥ 4.0)**
- **R packages:**
  - DESeq2
  - pheatmap
  - ggplot2
  - EnhancedVolcano (optional)
  - readr
  - dplyr
- **Shell utilities (bash, conda, etc.)** for environment setup

## Usage

### 1. Set up the environment

```bash
bash RNA_script.sh
```

This script installs dependencies and prepares your environment.

### 2. Run DESeq2 analysis

```r
source("DESeq2.r")
```

- Input: Raw count matrix and sample metadata
- Output: Normalized counts and a list of differentially expressed genes

### 3. Generate Heatmap

```r
source("heatmap.r")
```

- Input: DESeq2 results
- Output: Clustered heatmap of top differentially expressed genes

### 4. Create Volcano Plot

```r
source("volcano_plot.r")
```

- Input: DESeq2 results
- Output: Volcano plot highlighting significant genes

## Author

Assem Kadry Elsherif  
Assistant Lecturer, School of Biotechnology, Nile University  
Email: akadry@nu.edu.eg

## License

This project is licensed under the MIT License.
