# RNA-seq-spikein

This repository contains a customised data analysis pipeline that facilitates the incorporation of ERCC spike-in control for RNA-seq data analysis (For paired-end short read - Illumina data). 

### Table of contents 
1. Outline
2. Software installation
3. How to run - example dataset
  
     3.1 Quality check of raw fastq files - FastQC/MultiQC
  
     3.2 Read alignment, processing and post-alignment quality check - Bowtie2, Samtools and BamQC/MultiQC

     3.3 Coverage track generation and visualisation - deepTools and IGV
  
     3.4 Read summarisation - Subread package 'featureCounts'
  
     3.5 Differential expression analysis - R Bioconductor package 'DESeq2'
 
     3.6 Functional enrichment analysis - R Bioconductor package 'clusterProfiler'

     3.7 Gene level TPM estimation - Kallisto
  

### 1. Outline

This pipeline was run on an HPC (high-performance computing) system based on CentOS (Linux). R based analysis was carried out in RStudio/MacOS Catalina.
The 'Scripts' folder contains template bash scripts and example Snakemake workflows.

![alt text](https://github.com/tudumanne/RNA-seq-spikein/files/7829892/Picture.1.pdf)


2. Software installation 

The required software/command-line tools were installed via conda on Linux. 
(Miniconda https://docs.conda.io/en/latest/miniconda.html)


3. How to run - example dataset

The folder 'example dataset' contains 12 RNA-seq samples, 4 biological replicates per each stage (WT, PreM and Mal).
These files contain a subset of reads 1 million.  
  
3.1 Quality check of raw fastq files - FastQC/MultiQC
  
3.2 Read alignment, processing and post-alignment quality check - Bowtie2, Samtools and BamQC/MultiQC

3.3 Coverage track generation and visualisation - deepTools and IGV
  
3.4 Read summarisation - Subread package 'featureCounts'
  
3.5 Differential expression analysis - R Bioconductor package 'DESeq2'
 
3.6 Functional enrichment analysis - R Bioconductor package 'clusterProfiler'

3.7 Gene level TPM estimation - Kallisto
