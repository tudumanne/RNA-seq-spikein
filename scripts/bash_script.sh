#!/bin/bash

# usage - RNA-spikein analysis

#change the working directory
cd ~/rna-seq/experiment_01

#quality check of raw fastq bamfiles
fastqc -o fastqc --extract --dir fastqc --format fastq sample_*.fastq.gz
multiqc fastqc/

#read alignment and file processing
hisat2-build reference_genome.fa reference_genome+rDNA+ERCC
hisat2-align -q -x reference_genome+rDNA+ERCC -1 {sample}_R1.fastq -2 {sample}_R2.fastq --no-spliced-alignment --rna-strandness RF --time --no-unal | samtools view -bS > {sample}.bam
samtools sort {sample}.bam > {sample}_sorted.bam
samtools index {sample}_sorted.bam

#quality check of bamfiles
qualimap bamqc -bam file.bam -outdir qualimap_results -outformat pdf
multiqc fastqc/

#deeptools coverage track generation
bamCoverage -b {sample}_sorted.bam -o {sample}_coverage.bw --outFileFormat bigwig --binSize 10 --smoothLength 25 --extendReads --ignoreForNormalization chrX chrY BK000964.3 chrMT

#featurecounts read summarisation
featureCounts -p -a reference_genome.chr+ERCC.gtf -s 2 -t CDS -B -C -g gene_id -o counts_CDS_splice.txt *_sorted.bam

#kallisto tpm estimation
kallisto quant -i Mus_musculus.GRCm38.cdna.idx -o /sample {sample}_R1.fastq.gz {sample}_R2.fastq.gz

