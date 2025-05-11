
# Data Description
#This dataset contains paired-end RNA sequencing (RNA-seq) data from six biological samples, categorized into:
#
#Cancer Samples (3 biological replicates)
#    •   cancer_sample_1.read1.fastq.gz & read2.fastq.gz
#    •   cancer_sample_2.read1.fastq.gz & read2.fastq.gz
#    •   cancer_sample_3.read1.fastq.gz & read2.fastq.gz
#
#Each pair (read1 and read2) represents paired-end reads from a single cancer sample.
#
#Normal Samples (3 biological replicates)
#    •   normal_sample_1.read1.fastq.gz & read2.fastq.gz
#    •   normal_sample_2.read1.fastq.gz & read2.fastq.gz
#    •   normal_sample_3.read1.fastq.gz & read2.fastq.gz
#
#Each pair represents matched normal (non-cancerous) tissue from control samples.
#
#
#################################################################################################################


#### Set up directories
mkdir -p ~/workdir/{fqData,sample_data,trimmed,bwa_align/bwaIndex}


####  Install necessary software

mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh


~/miniconda3/bin/conda init bash

## restart the terminal
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda

## create the conda env if you have not already did
conda create -y --name ngs python=3.12
conda activate ngs
conda install -c bioconda fastqc multiqc bwa samtools trimmomatic subread
conda install conda-forge::gdown

#################################################################################################################

# Download Fastq files
cd ~/workdir/fqData
#gdown https://drive.google.com/uc?id=108ibdFuzypk_dNJa40pbnljL7zdXOkTB
#or 
#gdown 108ibdFuzypk_dNJa40pbnljL7zdXOkTB

tar -xzf RNA_seq_data.tar.gz

# Run FastQC on each fastq file
for f in *.fastq.gz; do
    fastqc -t 8 -f fastq -noextract $f
done

# Merge FastQC reports
multiqc -z -o . .

#################################################################################################################

# Error Trimming
cd ~/workdir/trimmed 

adap="$CONDA_PREFIX/share/trimmomatic-0.39-2/adapters"

# Perform error trimming on all raw data
for condition in cancer normal; do
  for i in 1 2 3; do
    read1=~/workdir/fqData/${condition}_sample_${i}.read1.fastq.gz
    read2=~/workdir/fqData/${condition}_sample_${i}.read2.fastq.gz

    trimmed_read1=~/workdir/trimmed/${condition}_sample_${i}.read1.trimmed.fastq.gz
    trimmed_read2=~/workdir/trimmed/${condition}_sample_${i}.read2.trimmed.fastq.gz

    # Run Trimmomatic to trim low-quality bases and adapters
    trimmomatic PE -phred33 $read1 $read2 \
      $trimmed_read1 ~/workdir/trimmed/${condition}_sample_${i}.read1.unpaired.fastq.gz \
      $trimmed_read2 ~/workdir/trimmed/${condition}_sample_${i}.read2.unpaired.fastq.gz \
      ILLUMINACLIP:$adap/TruSeq3-PE-2.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36
  done
done

#################################################################################################################

# Run Fastqc and multiqc on the trimmed data
for f in *trimmed.fastq.gz; do
    fastqc -t 8 -f fastq -noextract $f
done

# Merge FastQC reports
multiqc -z -o . .

#################################################################################################################

# Download reference genome and GTF file
cd ~/workdir/sample_data
#gdown 13_y4EMXQEI1weI2kGc7PaF4jIGKoCwvW
#gdown 1EeQLuJdKaNxb31GjHItazF4SzOtqHX-N

# Index the genome
bwa index -a bwtsw -p ~/workdir/bwa_align/bwaIndex/reference_genome_chr22 ~/workdir/sample_data/reference_genome_chr22.fa

bwaIndex=~/workdir/bwa_align/bwaIndex/reference_genome_chr22
GTF=~/workdir/sample_data/chr22.gtf
READS_DIR=~/workdir/fqData

# Create working directories
mkdir -p ~/workdir/diff_exp/{bams,sam_files}
cd ~/workdir/diff_exp


# Alignment Loop (using trimmed data)
for condition in cancer normal; do
  for i in 1 2 3; do
    read1=~/workdir/trimmed/${condition}_sample_${i}.read1.trimmed.fastq.gz
    read2=~/workdir/trimmed/${condition}_sample_${i}.read2.trimmed.fastq.gz
    sample=${condition}_${i}

    sam=~/workdir/diff_exp/sam_files/${sample}.sam
    bam=~/workdir/diff_exp/bams/${sample}.bam

    # Run BWA MEM for alignment
    bwa mem $bwaIndex $read1 $read2 > $sam

    # Convert SAM to sorted BAM
    samtools view -Sb $sam | \
    samtools sort -o $bam -

    # Index BAM
    samtools index $bam
  done
done

#################################################################################################################

# Quantification using featureCounts
featureCounts -p -a $GTF -g gene_name -o counts.txt  bams/cancer_*.bam  bams/normal_*.bam
# Simplify the file to keep only the count columns.
cat counts.txt | cut -f 1,7-12 > simple_counts.txt

# For differential expression, we will use DESeq R package and for visualization, we will use gplots package. 
conda install r
conda install r-gplots
conda install -c bioconda bioconductor-deseq2

mkdir -p ~/workdir/scripts && cd ~/workdir/scripts
wget https://raw.githubusercontent.com/assemkadry/NGS-R-Scripts/main/DESeq2.r # you can use curl -O instead of wget
wget https://raw.githubusercontent.com/assemkadry/NGS-R-Scripts/main/volcano_plot.r # you can use curl -O instead of wget
wget https://raw.githubusercontent.com/assemkadry/NGS-R-Scripts/main/heatmap.r # you can use curl -O instead of wget


## Differential expression by DESeq2
cd ~/workdir/diff_exp
cat simple_counts.txt | Rscript ~/workdir/scripts/DESeq2.r 3x3 > results_deseq2.tsv

#View only rows with padj < 0.05
cat results_deseq2.tsv | awk 'NR==1 || ($7 != "NA" && $7 < 0.05)' > filtered_results_deseq2.tsv

#Draw Volcan Plot
cat filtered_results_deseq2.tsv | Rscript ~/workdir/scripts/volcano_plot.r > volcano_plot.pdf
#Draw Heatmap
cat norm-matrix-deseq2.tsv | Rscript ~/workdir/scripts/heatmap.r > heatmap.pdf

