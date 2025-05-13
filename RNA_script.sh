#!/bin/bash
set -euo pipefail

# ===============================================================
# RNA-seq Paired-End Analysis Pipeline
# ===============================================================

# ---------------------------------------------------------------
# Data Descriptipn
# ---------------------------------------------------------------
# This dataset contains paired-end RNA sequencing (RNA-seq) data from six biological samples, categorized into:

# Cancer Samples (3 biological replicates)
#     •   cancer_sample_1.read1.fastq.gz & read2.fastq.gz
#     •   cancer_sample_2.read1.fastq.gz & read2.fastq.gz
#     •   cancer_sample_3.read1.fastq.gz & read2.fastq.gz

# Each pair (read1 and read2) represents paired-end reads from a single cancer sample.

# Normal Samples (3 biological replicates)
#     •   normal_sample_1.read1.fastq.gz & read2.fastq.gz
#     •   normal_sample_2.read1.fastq.gz & read2.fastq.gz
#     •   normal_sample_3.read1.fastq.gz & read2.fastq.gz

# Each pair represents matched normal (non-cancerous) tissue from control samples.

# ------------------------------
# Step 1: Set up directories
# ------------------------------
mkdir -p "$HOME"/workdir/{fqData,sample_data,trimmed,bwa_align/bwaIndex}

# ------------------------------
# Step 2: Install software with conda
# ------------------------------
mkdir -p "$HOME"/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O "$HOME"/miniconda3/miniconda.sh
bash "$HOME"/miniconda3/miniconda.sh -b -u -p "$HOME"/miniconda3
rm "$HOME"/miniconda3/miniconda.sh

"$HOME"/miniconda3/bin/conda init bash

# Restart the terminal
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda

# Create the conda env if you have not already did
conda create -y --name ngs python=3.12
conda activate ngs
conda install -c bioconda fastqc multiqc bwa samtools trimmomatic subread
conda install conda-forge::gdown

# ------------------------------
# Step 3: Download and extract FASTQ files
# ------------------------------
cd "$HOME"/workdir/fqData
#gdown 108ibdFuzypk_dNJa40pbnljL7zdXOkTB # Uncomment to download via gdown
tar -xzf RNA_seq_data.tar.gz # Extract FASTQ files

# ------------------------------
# Step 4: Run FastQC and MultiQC on raw FASTQ files
# ------------------------------
for f in *.fastq.gz; do
    fastqc -t 8 -f fastq -noextract "$f" # Run FastQC
done
multiqc -z -o . . # Aggregate FastQC reports

# ------------------------------
# Step 5: Adapter and quality trimming
# ------------------------------
cd "$HOME"/workdir/trimmed
adap="$CONDA_PREFIX/share/trimmomatic-0.39-2/adapters"

# Loop through samples and trim
for condition in cancer normal; do
  for i in 1 2 3; do
    read1="$HOME"/workdir/fqData/${condition}_sample_${i}.read1.fastq.gz
    read2="$HOME"/workdir/fqData/${condition}_sample_${i}.read2.fastq.gz

    trimmed_read1="$HOME"/workdir/trimmed/${condition}_sample_${i}.read1.trimmed.fastq.gz
    trimmed_read2="$HOME"/workdir/trimmed/${condition}_sample_${i}.read2.trimmed.fastq.gz

    # Run Trimmomatic to trim low-quality bases and adapters
    trimmomatic PE -phred33 -trimlog trimLogFile -summary statsSummaryFile $read1 $read2 \
      $trimmed_read1 ~/workdir/trimmed/${condition}_sample_${i}.read1.unpaired.fastq.gz \
      $trimmed_read2 ~/workdir/trimmed/${condition}_sample_${i}.read2.unpaired.fastq.gz \
      ILLUMINACLIP:$adap/TruSeq3-PE-2.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36  done
    done  
done

# ------------------------------
# Step 6: FastQC and MultiQC on trimmed reads
# ------------------------------
for f in *trimmed.fastq.gz; do
    fastqc -t 8 -f fastq -noextract "$f" # Run FastQC on trimmed reads
done
multiqc -z -o . . # Aggregate trimmed FastQC reports

# ------------------------------
# Step 7: Download reference genome and GTF
# ------------------------------
cd "$HOME"/workdir/sample_data
#gdown 13_y4EMXQEI1weI2kGc7PaF4jIGKoCwvW # Uncomment to use gdown
#gdown 1EeQLuJdKaNxb31GjHItazF4SzOtqHX-N # Uncomment to use gdown

bwa index -a bwtsw -p "$HOME"/workdir/bwa_align/bwaIndex/reference_genome_chr22 "$HOME"/workdir/sample_data/reference_genome_chr22.fa

bwaIndex="$HOME"/workdir/bwa_align/bwaIndex/reference_genome_chr22
GTF="$HOME"/workdir/sample_data/chr22.gtf

# ------------------------------
# Step 8: Align reads with BWA 
# ------------------------------
mkdir -p "$HOME"/workdir/diff_exp/{bams,sam_files}
cd "$HOME"/workdir/diff_exp

for condition in cancer normal; do
  for i in 1 2 3; do
    read1="$HOME"/workdir/trimmed/${condition}_sample_${i}.read1.trimmed.fastq.gz
    read2="$HOME"/workdir/trimmed/${condition}_sample_${i}.read2.trimmed.fastq.gz
    sample=${condition}_${i}

    sam="$HOME"/workdir/diff_exp/sam_files/${sample}.sam
    bam="$HOME"/workdir/diff_exp/bams/${sample}.bam

    # Align with BWA MEM
    bwa mem "$bwaIndex" "$read1" "$read2" > "$sam"

    # Convert SAM to sorted BAM
    samtools view -Sb "$sam" | samtools sort -o "$bam" - 
    # Index BAM
    samtools index "$bam"
  done
done

# ------------------------------
# Step 9: Gene quantification with featureCounts
# ------------------------------
featureCounts -p -a "$GTF" -g gene_name -o counts.txt bams/cancer_*.bam bams/normal_*.bam
# Simplify the file to keep only the count columns.
cut -f 1,7-12 counts.txt > simple_counts.txt

# ------------------------------
# Step 10: Differential expression with DESeq2
# ------------------------------
# For differential expression, we will use DESeq2 R package and for visualization, we will use gplots package. 
conda install r
conda install r-gplots
conda install -c bioconda bioconductor-deseq2

mkdir -p "$HOME"/workdir/scripts && cd "$HOME"/workdir/scripts
wget https://raw.githubusercontent.com/assemkadry/NGS-R-Scripts/main/DESeq2.r
wget https://raw.githubusercontent.com/assemkadry/NGS-R-Scripts/main/volcano_plot.r
wget https://raw.githubusercontent.com/assemkadry/NGS-R-Scripts/main/heatmap.r

cd "$HOME"/workdir/diff_exp
# Run DESeq2
cat simple_counts.txt | Rscript "$HOME"/workdir/scripts/DESeq2.r 3x3 > results_deseq2.tsv

# View only rows with padj < 0.05
awk 'NR==1 || ($7 != "NA" && $7 < 0.05)' results_deseq2.tsv > filtered_results_deseq2.tsv

# Draw Volcan Plot
cat filtered_results_deseq2.tsv | Rscript "$HOME"/workdir/scripts/volcano_plot.r > volcano_plot.pdf
# Draw Heatmap
cat norm-matrix-deseq2.tsv | Rscript "$HOME"/workdir/scripts/heatmap.r > heatmap.pdf
