#functional enrichment analysis of RNA-seq data - identified as differentially expressed
#R Bioconductor package 'clusterProfiler'

#Reference
#https://yulab-smu.top/biomedical-knowledge-mining-book/index.html

#load libraries 
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggnewscale)
library(DOSE)
library(ReactomePA)

#read in gene list 
data1 <- read.csv("wt_mal_ercc_results_padj0.05.csv", header = TRUE)

#set the desired organism
organism = org.Mm.eg.db

#log2 fold change as a vector
original_gene_list1 <- data1$log2FoldChange

#name the vector
names(original_gene_list1) <- data1$X

#omit any NA values 
gene_list1<-na.omit(original_gene_list1)

#sort the list in decreasing order (required for clusterProfiler)
gene_list1 = sort(gene_list1, decreasing = TRUE)

#gene set enrichment analysis (gsea) of gene ontology (GO) categories
gse1 <- gseGO(geneList=gene_list1, 
              ont ="BP", 
              keyType = "ENSEMBL", 
              minGSSize = 10, 
              maxGSSize = 500, 
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = organism, 
              pAdjustMethod = "none",
              eps = 0)

#dot plot of GO terms
dotplot(gse1, showCategory=15) + facet_grid(.~.sign)

#plot distribution of core enriched genes for GSEA enriched categories
ridgeplot(gse1, showCategory = 15) + labs(x = "enrichment distribution")+xlim(-10,10)

#plot running score and pre-ranked list of genes for a specific enriched category
gseaplot2(gse1, geneSetID = 1, title = gse1$Description[1])
