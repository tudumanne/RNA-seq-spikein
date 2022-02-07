# RNA-seq-spikein

This repository contains a customised data analysis pipeline that facilitates the incorporation of ERCC spike-in control for RNA-seq data analysis (For paired-end short read - Illumina data). 

### Table of contents 
1. [Overview of the pipeline](#overview-of-the-pipeline)
2. [Software installation](#software-installation)
3. [How to run an example dataset](#how-to-run-an-example-dataset)
  
     3.1 Quality check of raw fastq files - FastQC/MultiQC
  
     3.2 Read alignment, processing and post-alignment quality check - HISAT2, Samtools and BamQC/MultiQC

     3.3 Coverage track generation and visualisation - deepTools and IGV
  
     3.4 Read summarisation - Subread package 'featureCounts'
  
     3.5 Differential expression analysis - R Bioconductor package 'DESeq2'
 
     3.6 Functional enrichment analysis - R Bioconductor package 'clusterProfiler'

     3.7 Gene level TPM estimation - Kallisto
  

### Overview of the pipeline

This pipeline was run on an HPC (high-performance computing) system based on CentOS (Linux). R based analysis was carried out in RStudio/MacOS Catalina.
The 'Scripts' folder contains template bash scripts for the analysis.

![alt text](https://github.com/tudumanne/RNA-seq-spikein/files/7829892/Picture.1.pdf)


### Software installation 

The required software/command-line tools were installed via conda on Linux. 
(Miniconda https://docs.conda.io/en/latest/miniconda.html)

```console
conda env create -n chip-seq -f environment.yml
```

### How to run an example dataset

The folder 'example dataset' contains 12 RNA-seq samples, 4 biological replicates per each stage (WT, PreM and Mal).
These files contain a subset of reads 1 million.  
  
3.1 Quality check of raw fastq files - FastQC/MultiQC

```console
fastqc -o fastqc --extract --dir fastqc --format fastq sample_*.fastq.gz
multiqc fastqc/
```
  
3.2 Read alignment, processing and post-alignment quality check - HISAT2, Samtools and BamQC/MultiQC

```console
hisat2-build reference_genome.fa reference_genome+rDNA+ERCC
hisat2-align -q -x reference_genome+rDNA+ERCC -1 {sample}_R1.fastq -2 {sample}_R2.fastq --no-spliced-alignment --rna-strandness RF --time --no-unal | samtools view -bS > {sample}.bam
samtools sort {sample}.bam > {sample}_sorted.bam
samtools index {sample}_sorted.bam
```

Quality check of BAM files

```console
qualimap bamqc -bam file.bam -outdir qualimap_results -outformat pdf
multiqc fastqc/
``` 

3.3 Coverage track generation and visualisation - deepTools and IGV

Generate a coverage track

```console
bamCoverage -b {sample}_sorted.bam -o {sample}_coverage.bw --outFileFormat bigwig --binSize 10 --smoothLength 25 --extendReads --ignoreForNormalization chrX chrY BK000964.3 chrMT 
```
  
3.4 Read summarisation - Subread package 'featureCounts'

```console
featureCounts -p -a reference_genome.chr+ERCC.gtf -s 2 -t CDS -B -C -g gene_id -o counts_CDS_splice.txt *_sorted.bam
```
  
3.5 Differential expression analysis - R Bioconductor package 'DESeq2'

-R script-
 
3.6 Functional enrichment analysis - R Bioconductor package 'clusterProfiler'

-R script-

3.7 Gene level TPM estimation - Kallisto

```console
kallisto quant -i Mus_musculus.GRCm38.cdna.idx -o /sample {sample}_R1.fastq.gz {sample}_R2.fastq.gz
```

