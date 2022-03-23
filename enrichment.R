if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msigdf")
library(DESeq2)
rm(list=ls())
gene_count <- read.csv("counts.csv",stringsAsFactors = F)
library(stringr)
gene_count_mod <- rep(NA,nrow(gene_count))
gene_count$gene_id <-gene_count_mod
rownames(gene_count) <- gene_count[,1]
gene_count <- gene_count[,-1]
gene_count_CON <- gene_count[,c(1,2)]
gene_count_Yoda <- gene_count[,c(3,4)]
gene_count_use <- gene_count[,c(1,2,3,4)]
write.csv(gene_count_control,file="control")
write.csv(gene_count_OE,file="Yoda")
condition<- factor(c(rep("CON"),rep("CON"),rep("Yoda"),rep("Yoda")),levels=c("CON","Yoda"))
condition
colData <- data.frame(row.names = colnames(gene_count_use),condition)
colData
dds_Yoda <- DESeqDataSetFromMatrix(gene_count_use,colData,design = ~condition)
dds_Yoda <- DESeq(dds_Yoda)
dds_Yoda
res_Yoda = results(dds_Yoda, contrast = c("condition","CON","Yoda"))
res_Yoda = res_Yoda[order(res_Yoda$pvalue),]
summary(res_Yoda)
write.csv(res_Yoda,file="Def_Expr_Genes.csv")

library(ggplot2)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
BiocManager::install("clusterProfiler",force = TRUE)
library(clusterProfiler)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
BiocManager::install("org.Mm.eg.db")
BiocManager::install("DOSE")
BiocManager::install("topGO")
library(org.Mm.eg.db)
library(DOSE)
library(topGO)
keytypes(org.Mm.eg.db) 

T4_0[which(T4_0$T4vsT0_padj %in% NA),'sig'] <- 'no diff'
T4_0[which(T4_0$T4vsT0_log2FoldChange >= 1 & T4_0$T4vsT0_padj < 0.05),'sig'] <- 'up'
T4_0[which(T4_0$T4vsT0_log2FoldChange <= -1 & T4_0$T4vsT0_padj < 0.05),'sig'] <- 'down'
T4_0[which(abs(T4_0$T4vsT0_log2FoldChange) < 1 | T4_0$T4vsT0_padj >= 0.05),'sig'] <- 'no diff'
data_up <- subset(T4_0, sig=='up')
data_down <- subset(T4_0, sig=="down")
data <- rbind(data_up, data_down)
all_gene <- c()
for (i in data$gene_id){
  print(i)
  all_gene <- c(all_gene,i)
  
}
all_gene_list = bitr(all_gene, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID",'UNIPROT'), OrgDb="org.Mm.eg.db")
ego_all <- enrichGO(gene = all_gene_list$ENTREZID, OrgDb = org.Mm.eg.db,ont = "ALL")
head(ego_all)
ego_all_1 <- setReadable(ego_all, OrgDb = org.Mm.eg.db)
KEGG<-enrichKEGG(gene=all_gene_list$ENTREZID,organism = 'mouse',pvalueCutoff = 0.9,qvalueCutoff = 0.9)
kk=DOSE::setReadable(KEGG,OrgDb = 'org.Mm.eg.db',keyType = 'ENTREZID')
dotplot(kk,title="T4_All_DEGs_KEGG")

