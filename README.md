# 🧬 RNA-Seq Differential Expression Pipeline

> A reproducible RNA-seq analysis pipeline for educational and research use

---

This repository provides a structured pipeline for analyzing paired-end RNA-Seq data. It covers preprocessing, alignment, quantification, and differential expression analysis using open-source tools.  
The pipeline was developed as part of the **Next Generation Sequencing (NGS)** course at the **School of Biotechnology, Nile University, Egypt**.

---

## 📑 Table of Contents
- [Educational Context](#educational-context)
- [Dataset Overview](#dataset-overview)
- [Pipeline Steps](#pipeline-steps)
- [Tools Used](#️tools-used)
- [Folder Structure](#folder-structure)
- [How to Run](#how-to-run)
- [Output](#output)
- [System Requirements](#system-requirements)
- [Example Output Files](#example-output-files)
- [Instructor](#instructor)
- [License](#license)

---

## 📚 Educational Context

This pipeline was designed for undergraduate students enrolled in the NGS course.  
It provides hands-on experience with standard tools and workflows in transcriptomics.

---

## 📁 Dataset Overview

**Samples:**
- 🧪 **Cancer (3 replicates)**  
  `cancer_sample_1.read1.fastq.gz` & `read2.fastq.gz`  
  `cancer_sample_2.read1.fastq.gz` & `read2.fastq.gz`  
  `cancer_sample_3.read1.fastq.gz` & `read2.fastq.gz`  

- 🧬 **Normal (3 replicates)**  
  `normal_sample_1.read1.fastq.gz` & `read2.fastq.gz`  
  `normal_sample_2.read1.fastq.gz` & `read2.fastq.gz`  
  `normal_sample_3.read1.fastq.gz` & `read2.fastq.gz`

---

## 🧪 Pipeline Steps

1. **Directory Setup**  
   Create folders for raw data, trimmed reads, reference genome, and alignments.

2. **Tool Installation**  
   Install Miniconda and required packages using Bioconda.

3. **Quality Control**  
   Run FastQC and summarize results with MultiQC.

4. **Trimming**  
   Trim low-quality bases and adapters using Trimmomatic.

5. **Alignment**  
   - Index the reference genome using BWA  
   - Align reads using BWA MEM  
   - Convert and sort SAM to BAM using SAMtools

6. **Quantification**  
   Use Subread’s `featureCounts` to count mapped reads per gene.

7. **Differential Expression Analysis**  
   Use DESeq2 in R to identify differentially expressed genes.

8. **Visualization**  
   Generate volcano plots and heatmaps using R scripts.

---

## 🛠️ Tools Used

| Step                | Tool                      | Description                                  |
|---------------------|---------------------------|----------------------------------------------|
| QC                  | `FastQC`, `MultiQC`        | Raw read quality assessment                  |
| Trimming            | `Trimmomatic`              | Adapter and base quality trimming            |
| Alignment           | `BWA`, `SAMtools`          | Read alignment and sorting                   |
| Quantification      | `Subread (featureCounts)`  | Gene-level read counting                     |
| Analysis            | `R`, `DESeq2`              | DEG analysis and statistical testing         |
| Visualization       | `ggplot2`, `pheatmap`      | Volcano and heatmap plotting in R            |

---

## 📁 Folder Structure

```
~/workdir/
├── fqData/          # Raw FASTQ files
├── trimmed/         # Cleaned reads
├── sample_data/     # Reference genome & GTF
├── bwa_align/
│   ├── bwaIndex/    # BWA index files
│   └── alignments/  # Sorted BAM files
```

---

## ▶️ How to Run

```bash
bash RNA_script.sh
```

✅ Ensure the reference genome (FASTA) and annotation (GTF) are placed in `sample_data/`.

---

## 🧾 Output

- 📑 Quality Control Reports (FastQC/MultiQC)
- ✂️ Trimmed FASTQ files
- 🧬 Aligned BAM files
- 🧮 Counts matrix (from featureCounts)
- 📊 DEG results (DESeq2)
- 🖼 Volcano plot and heatmap

---

## 💻 System Requirements

- Unix/Linux system
- 8GB+ RAM
- Internet connection (for installing dependencies)
- Disk space: ~20 GB recommended

---

## 📂 Example Output Files

- `counts_matrix.txt`
- `deseq2_results.csv`
- `volcano_plot.pdf`
- `heatmap.pdf`

---

## 👨‍🏫 Instructor

**Assem Kadry Elsherif**  
Assistant Lecturer, School of Biotechnology  
Nile University, Egypt

📧 akadry@nu.edu.eg  
🌍 GitHub: [assem-kadry](https://github.com/assem-kadry)

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
