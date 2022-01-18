#Differential expression analysis of RNA-seq data
#R Bioconductor package 'deseq2'

#Reference
#https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt
#https://www.biostars.org/p/166556/
#https://www.biostars.org/p/166556/

#input file - featurecounts output file
#separate .txt file containing ERCC spike-in counts

#load libraries
library(ggpubr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(vsn)
library(pheatmap)
library("RColorBrewer")
library(GenomicFeatures)
library(org.Mm.eg.db)
library("AnnotationDbi")
library(EnsDb.Mmusculus.v79)

#loading ercc spikein counts and renaming columns
countdata_ercc <- read.table("counts_exon_splice_ercc.txt", header=TRUE, row.names=1)
countdata_ercc <- countdata_ercc[ ,6:ncol(countdata_ercc)]
colnames(countdata_ercc) <- gsub("\\.[sb]am$", "", colnames(countdata_ercc))
colnames(countdata_ercc) <- gsub("\\_sorted$", "", colnames(countdata_ercc))
#count data as a matrix
countdata_ercc <- as.matrix(countdata_ercc)
head(countdata_ercc)

#define the different conditions (stages) and biological replicates
(stage <- factor(c(rep("wt", 4), rep("prem", 4), rep("mal", 4))))
#as a data frame
(coldata_ercc <- data.frame(row.names=colnames(countdata_ercc), stage))

#store count data in a deseqdataset object
dds_ercc <- DESeqDataSetFromMatrix(countData=countdata_ercc, colData=coldata_ercc, design=~stage)
dds_ercc
#estimate size factors for each sample based on ercc spike-in counts
dds_ercc <- estimateSizeFactors(dds_ercc)
dds_ercc
#assign size factors to a variable 
nf <- sizeFactors(dds_ercc)
sizeFactors(dds_ercc)

#loading rest of the data and renaming the columns
countdata <- read.table("counts_exon_splice_noercc.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("\\_sorted$", "", colnames(countdata))
#count data as a matrix
countdata <- as.matrix(countdata)
head(countdata)

(coldata <- data.frame(row.names=colnames(countdata), stage))

#store count data in a deseqdataset object
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~stage)
dds
#assign size factors calculated based on ercc (instead of calculating it based on rest of the data)
sizeFactors(dds) <- nf

#filtering rows with counts less that 10
keep <- rowSums((counts(dds)) >= 10) >= 4
dds <- dds[keep,]

#run differential expression analysis 
dds <- DESeq(dds)
res <- results(dds)
#view results
res

#write normalised count data to a .csv file
data_table = as.data.frame(counts(dds, normalized=TRUE))
write.csv(as.data.frame(data_table), 
          file="countdata_all.csv")

#define contrasts WT-PreM, PreM-Mal, WT-Mal 
#e.g. WT-Mal
res.W.M <- results(dds, contrast=c("stage","mal", "wt"))
#view genes that are significantly differentially expressed (p.adj<0.05)
table(res.W.M$padj<0.05)
res.W.M <- res.W.M[order(res.W.M$padj), ]
#histogram of p-values
hist( res.W.M$pvalue, breaks=20, col="grey" )


#plot the normalised counts of a selected gene
d <- plotCounts(dds, gene="ENSMUSG00000022346", intgroup="stage", returnData=TRUE)
d$stage <- factor(d$stage, levels = c("wt", "prem", "mal"))

ggplot(d, aes(x = stage, y = count, color=stage)) + 
  geom_point(position=position_jitter(w = 0.1, h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("Myc - normalized counts") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(colour = "black"), plot.title=element_text(face="bold"), plot.subtitle=element_text(face="bold"), plot.caption =element_text(face="bold"),
        axis.title.x=element_text(face="bold", size=12, colour = "black"), axis.text.x = element_text(face = "bold", size=12, colour = "black"),
        axis.title.y=element_text(face="bold", size=12, colour = "black"), axis.text.y = element_text(face = "bold", size=12, colour = "black"), 
        legend.text = element_text(face = "bold", size=12, colour = "black"), legend.title = element_text(face = "bold", size=12, colour = "black") )

#heatmap of sample-to-sample distance
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rownames(coldata))
colnames(sampleDistMatrix) <- paste(rownames(coldata))
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#plot PCA
pcaData <- plotPCA(rld, intgroup=c("stage"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=stage)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  geom_label_repel(aes(label = rownames(pcaData), size = 1)) +
  theme_bw()+
  theme(text = element_text(colour = "black"), plot.title=element_text(face="bold"), plot.subtitle=element_text(face="bold"), plot.caption =element_text(face="bold"),
        axis.title.x=element_text(face="bold", size=12, colour = "black"), axis.text.x = element_text(face = "bold", size=12, colour = "black"),
        axis.title.y=element_text(face="bold", size=12, colour = "black"), axis.text.y = element_text(face = "bold", size=12, colour = "black"), 
        legend.text = element_text(face = "bold", size=12, colour = "black"), legend.title = element_text(face = "bold", size=12, colour = "black") )

#plot dispersion of normalised data
plotDispEsts(dds)


#create a volcano plot with 'EnhancedVolcano' package
library(EnhancedVolcano)

EnhancedVolcano(res.W.M,
                lab = rownames(res.W.M),
                x = 'log2FoldChange',
                y = 'pvalue',
                # selectLab = c("ENSMUSG00000022346"), 
                labSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-10,15))


#export DEG list as a .csv file
write.csv(as.data.frame(res.W.M), 
          file="wt_mal_ercc_results.csv")

#only significant DEGs (padj<0.05)
res.WM.sig <- subset(res.W.M, padj < 0.05)
write.csv(as.data.frame(res.WM.sig), 
          file="wt_mal_ercc_results_padj0.05.csv")




