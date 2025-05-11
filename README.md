# 🧬 RNA-Seq Differential Expression Pipeline

> 🔬 *End-to-end workflow for preprocessing, aligning, counting, and analyzing paired-end RNA-seq data.*

---

## 📚 Educational Context

This pipeline was designed as a teaching resource for the **Next Generation Sequencing (NGS)** course  
at the **School of Biotechnology, Nile University, Egypt**.  
It guides undergraduate students through a practical RNA-seq data analysis workflow using real tools and datasets.

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

<pre><code>```mermaid
graph TD
  A[Start] --> B[Directory Setup]
  B --> C[Tool Installation (via Bioconda)]
  C --> D[Quality Control (FastQC + MultiQC)]
  D --> E[Trimming (Trimmomatic)]
  E --> F[Alignment (BWA + SAMtools)]
  F --> G[Counting (Subread/featureCounts)]
  G --> H[DESeq2 Analysis + Visualization]
  H --> I[Output Files]
```</code></pre>


---

## 🛠️ Tools Used

| Step                | Tool                      | Description                                  |
|---------------------|---------------------------|----------------------------------------------|
| QC                  | `FastQC`, `MultiQC`        | Raw read quality assessment                  |
| Trimming            | `Trimmomatic`              | Adapter and low-quality base removal         |
| Alignment           | `BWA`, `SAMtools`          | Read alignment and BAM sorting               |
| Read Quantification | `Subread (featureCounts)`  | Count reads mapped to genes                  |
| Analysis & Plots    | `R + DESeq2`               | DEG analysis, volcano & heatmap generation   |

---

## 📁 Folder Structure

```
~/workdir/
├── fqData/          # Raw FASTQ files
├── trimmed/         # Cleaned reads after Trimmomatic
├── sample_data/     # Reference genome & annotation
├── bwa_align/
│   ├── bwaIndex/    # Indexed reference files
│   └── alignments/  # Aligned BAM/SAM files
```

---

## ▶️ How to Run

```bash
bash RNA_script.sh
```

✅ Make sure your reference genome (FASTA) and annotation (GTF) are in `sample_data/`.

---

## 🧾 Output

- 📑 **QC Reports**: HTML summaries from FastQC and MultiQC  
- 📎 **Trimmed Reads**: Clean FASTQ files  
- 📂 **Alignments**: Sorted BAM files  
- 📊 **Counts Table**: Gene expression matrix (from `featureCounts`)  
- 📈 **Plots**: Volcano plot and heatmap (from R scripts)

---

## 💻 System Requirements

- 🐧 Unix/Linux OS
- 🧠 8GB+ RAM
- 🌐 Internet access
- 📦 Conda (installed via Miniconda in the script)
- 📊 R + Bioconductor packages

---

## 👨‍🏫 Instructor

Developed and maintained by  
**Assem Kadry Elsherif**  
Assistant Lecturer, School of Biotechnology  
Nile University, Egypt

📧 akadry@nu.edu.eg  
🌍 GitHub: [assem-kadry](https://github.com/assem-kadry)
