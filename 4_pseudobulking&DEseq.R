#- Project 'E:/seurat/DESeq2' loaded. [renv 1.1.4]
library(apeglm)
library(DESeq2)
library(Seurat)
library(presto)
library(org.Dr.eg.db)
library(ggplot2)

#=============================loading===========================================
Hb<-readRDS("C:/Users/admin/Desktop/tidy up code/Hb_allgenotype_withintegration.rds")


#=========================DEseq for whole habenula==============================
#create DEseq object
peudo<-AggregateExpression(Hb, group.by = c("orig.ident"),return.seurat = FALSE)
print(colnames(peudo[["RNA"]]))

sample_info <- data.frame(
  row.names = colnames(counts),
  genotype = c("Mut", "Mut", "WT", "WT")
)
sample_info$genotype <- factor(sample_info$genotype)
sample_info$genotype <- relevel(sample_info$genotype, ref = "WT")

dds <- DESeqDataSetFromMatrix(
  countData = peudo[["RNA"]],
  colData = sample_info,
  design = ~ genotype
)

#filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds<- dds[keep,]
#Running DEseq algorithm
dds<- DESeq(dds)
resultsNames(dds)
res <-results(dds, name = "genotype_Mut_vs_WT")
#save results
res_df <- as.data.frame(res[order(res$log2FoldChange),])  #order the result by pvalue (optional) or log2FoldChange
res_df$gene<-rownames(res_df)
#write.csv(res_df,file = 'C:/Users/admin/Desktop/Hb_Deseq.csv')  ##

#=============================DEseq for each cluster across genotype============================
cluster<-subset(Hb, subset=cluster=="4") #!!!change here for different clusters: 
#Hb02(1) Hb05(9) Hb06(6) Hb07(7) Hb10(11) Hb15(3)
#Hb04(5) Hb08(4) Hb13(8)

#create DEseq object
peudo<-AggregateExpression(cluster, group.by = c("orig.ident"),return.seurat = FALSE)
print(colnames(peudo[["RNA"]]))

sample_info <- data.frame(
  row.names = colnames(counts),
  genotype = c("Mut", "Mut", "WT", "WT")
)
sample_info$genotype <- factor(sample_info$genotype)
sample_info$genotype <- relevel(sample_info$genotype, ref = "WT")

dds <- DESeqDataSetFromMatrix(
  countData = peudo[["RNA"]],
  colData = sample_info,
  design = ~ genotype
)

#filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds<- dds[keep,]
#Running DEseq algorithm
dds<- DESeq(dds)
resultsNames(dds)
res <-results(dds, name = "genotype_Mut_vs_WT")
#save results
res_df <- as.data.frame(res[order(res$log2FoldChange),])  #order the result by pvalue (optional) or log2FoldChange
res_df$gene<-rownames(res_df)
#write.csv(res_df,file = 'C:/Users/admin/Desktop/Hb08_Deseq.csv')  ##


#=======================DEseq for each cluster in WT============================
set.seed(123)
WT$target <- ifelse(WT$cluster == 12, "target", "other") #!!!change here for different clusters: 
#Hb01(0) Hb02(1) Hb05(9) Hb06(6) Hb07(7) Hb09(2) Hb10(11) Hb15(3)
#Hb04(5) Hb08(4) Hb13(8) Hb03(10) Hb11 (12)
WT$target<-paste0(WT$orig.ident,WT$target)

peudo<-AggregateExpression(WT, group.by = c("target"),return.seurat = FALSE)
print(colnames(peudo[["RNA"]]))

#create DEseq object
sample_info <- data.frame(
  row.names = colnames(counts),
  genotype = c("other", "target", "other", "target")
)
sample_info$genotype <- factor(sample_info$genotype)
sample_info$genotype <- relevel(sample_info$genotype, ref = "other")

dds <- DESeqDataSetFromMatrix(
  countData = peudo[["RNA"]],
  colData = sample_info,
  design = ~ genotype
)

#filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds<- dds[keep,]
#Running DEseq algorithm
dds<- DESeq(dds)
resultsNames(dds)
res <-results(dds, name = "genotype_target_vs_other")
#save results
res_df <- as.data.frame(res[order(res$log2FoldChange),])  #order the result by pvalue (optional) or log2FoldChange
res_df$gene<-rownames(res_df)
#write.csv(res_df,file = 'C:/Users/admin/Desktop/Hb11_WT_Deseq.csv')  ##













