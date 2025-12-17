#- Project 'E:/seurat/seurat_basics' loaded. [renv 1.1.4]
#[Workspace loaded from E:/seurat/seurat_basics/.RData]
library(Seurat)
library(presto)
library(svglite)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(SingleCellExperiment)
library(scater)
library(scran)
library(miQC)
library(DoubletFinder)


#==============================reading raw count matrix=======================================
WT1<-readRDS("C:/Users/admin/Desktop/seurat/2024-10-13 100K scseq/u775Tg-WT1-100K_Seurat.rds")
WT1$cell<-paste0("WT1_",rownames(WT1@meta.data))
WT1<-RenameCells(WT1, new.names = WT1$cell)
WT2<-readRDS("C:/Users/admin/Desktop/seurat/2024-10-13 100K scseq/u775Tg-WT2-100K_Seurat.rds")
WT2$cell<-paste0("WT2_",rownames(WT2@meta.data))
WT2<-RenameCells(WT2, new.names = WT2$cell)
Mut3<-readRDS("C:/Users/admin/Desktop/seurat/2024-10-13 100K scseq/u775Tg-Mut3-100K_Seurat.rds")
Mut3$cell<-paste0("Mut3_",rownames(Mut3@meta.data))
Mut3<-RenameCells(Mut3, new.names = Mut3$cell)
Mut4<-readRDS("C:/Users/admin/Desktop/seurat/2024-10-13 100K scseq/u775Tg-Mut4-100K_Seurat.rds")
Mut4$cell<-paste0("Mut4_",rownames(Mut4@meta.data))
Mut4<-RenameCells(Mut4, new.names = Mut4$cell)


#=================crude clustering and detection of habenular marker genes=============================
seurat_obj <-merge(x = WT1, y = c(WT2, Mut3, Mut4)) #or seperate libraries for appendix
seurat_obj$genotype<-ifelse(grepl("Mut", seurat_obj$orig.ident),"mutant","WT")

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.005)  #change here for different clustering resolution
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

fulldata <- seurat_obj

#ploting marker genes
#p1<-FeaturePlot(seurat_obj,c("pou4f1","irx1b"),reduction = 'umap') #habenula marker genes
#p2<-FeaturePlot(seurat_obj,c("gng8","lmo3"),reduction = 'umap')  #dHb and vHb marker genes
#p3<-FeaturePlot(seurat_obj,c("epcam","pvalb5"),reduction = 'umap') #olfactory neuron marker genes
#p4<-FeaturePlot(seurat_obj,c("elavl3","elavl4"),reduction = 'umap')  #nerve system marker genes
#p5a<-FeaturePlot(seurat_obj,c("EGFP"),reduction = 'umap')

#p5b<-DimPlot(seurat_obj, reduction = "pca", label = FALSE, group.by = "seurat_clusters") #ploting clusters 
#p5<-grid.arrange(p5a,p5b,ncol = 2)


#===============================filtering poor quality cells=======================================
#########method1: filtering poor quality cells on WT-mutant merged dataset using micQC 
seurat_obj <- as.SingleCellExperiment(merge(x = WT1, y = c(WT2, Mut3, Mut4)))
mt_genes<-rownames(seurat_obj)[grepl("^mt-", rownames(seurat_obj))]
seurat_obj<-addPerCellQC(seurat_obj,subsets = list(mito = mt_genes))
model<-mixtureModel(seurat_obj)
#plotModel(seurat_obj,model)
#plotFiltering(seurat_obj,model)
seurat_obj<-filterCells(seurat_obj,model)
seurat_obj<-as.Seurat(seurat_obj)

idx<-match(rownames(fulldata@meta.data),rownames(seurat_obj@meta.data))
fulldata@meta.data$miQC_mergeclust<-seurat_obj@meta.data[idx,"cell"]


############method2: filtering poor quality cells on seperate dataset using scater-MAD   
seurat_obj <- WT1
seurat_obj <- as.SingleCellExperiment(seurat_obj)
mt_genes<-rownames(seurat_obj)[grepl("^mt-", rownames(seurat_obj))]
seurat_obj<-addPerCellQC(seurat_obj,subsets = list(mito = mt_genes))
lib_out <- isOutlier(seurat_obj$sum, nmads = 3, type = "both")
feat_out <- isOutlier(seurat_obj$detected, nmads = 3, type = "both")
discard <- lib_out|feat_out
seurat_obj <- seurat_obj[, !discard]
seurat_obj <- as.Seurat(seurat_obj)
seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj,pattern = "^mt-")
seurat_obj<-subset(seurat_obj, subset=percent.mt <5)
#VlnPlot(seurat_obj,features = c("percent.mt","nFeature_RNA","nCount_RNA"))
WT1_filt<-seurat_obj

seurat_obj <- WT2  #run following lines 4 times (WT1, WT2, Mut3, Mut4) for each library!!!
seurat_obj <- as.SingleCellExperiment(seurat_obj)
mt_genes<-rownames(seurat_obj)[grepl("^mt-", rownames(seurat_obj))]
seurat_obj<-addPerCellQC(seurat_obj,subsets = list(mito = mt_genes))
lib_out <- isOutlier(seurat_obj$sum, nmads = 3, type = "both")
feat_out <- isOutlier(seurat_obj$detected, nmads = 3, type = "both")
discard <- lib_out|feat_out
seurat_obj <- seurat_obj[, !discard]
seurat_obj <- as.Seurat(seurat_obj)
seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj,pattern = "^mt-")
seurat_obj<-subset(seurat_obj, subset=percent.mt <5)
#VlnPlot(seurat_obj,features = c("percent.mt","nFeature_RNA","nCount_RNA"))
WT2_filt<-seurat_obj

seurat_obj <- Mut3  #run following lines 4 times (WT1, WT2, Mut3, Mut4) for each library!!!
seurat_obj <- as.SingleCellExperiment(seurat_obj)
mt_genes<-rownames(seurat_obj)[grepl("^mt-", rownames(seurat_obj))]
seurat_obj<-addPerCellQC(seurat_obj,subsets = list(mito = mt_genes))
lib_out <- isOutlier(seurat_obj$sum, nmads = 3, type = "both")
feat_out <- isOutlier(seurat_obj$detected, nmads = 3, type = "both")
discard <- lib_out|feat_out
seurat_obj <- seurat_obj[, !discard]
seurat_obj <- as.Seurat(seurat_obj)
seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj,pattern = "^mt-")
seurat_obj<-subset(seurat_obj, subset=percent.mt <5)
#VlnPlot(seurat_obj,features = c("percent.mt","nFeature_RNA","nCount_RNA"))
Mut3_filt<-seurat_obj

seurat_obj <- Mut4  #run following lines 4 times (WT1, WT2, Mut3, Mut4) for each library!!!
seurat_obj <- as.SingleCellExperiment(seurat_obj)
mt_genes<-rownames(seurat_obj)[grepl("^mt-", rownames(seurat_obj))]
seurat_obj<-addPerCellQC(seurat_obj,subsets = list(mito = mt_genes))
lib_out <- isOutlier(seurat_obj$sum, nmads = 3, type = "both")
feat_out <- isOutlier(seurat_obj$detected, nmads = 3, type = "both")
discard <- lib_out|feat_out
seurat_obj <- seurat_obj[, !discard]
seurat_obj <- as.Seurat(seurat_obj)
seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj,pattern = "^mt-")
seurat_obj<-subset(seurat_obj, subset=percent.mt <5)
#VlnPlot(seurat_obj,features = c("percent.mt","nFeature_RNA","nCount_RNA"))
Mut4_filt<-seurat_obj 

seurat_obj<-merge(x=WT1_filt,y=c(WT2_filt,Mut3_filt,Mut4_filt))
idx<-match(rownames(fulldata@meta.data),rownames(seurat_obj@meta.data))
fulldata@meta.data$MAD_seperateclust<-seurat_obj@meta.data[idx,"cell"]


#====================================doublet detection==============================
################method1:doublets detection without filtering
seurat_obj <-WT1 
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
#seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
#seurat_obj <- FindClusters(seurat_obj, resolution = 0.01)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
#optimization of pK. optimized pK value is applied to doubleFinder()
#sweep.res.list_seu <- paramSweep(seurat_obj, PCs = 1:30, sct = FALSE)
#sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
#bcmvn_seu <- find.pK(sweep.stats_seu)
#saveRDS(sweep.res.list_seu,file = "C:/Users/admin/Desktop/sweep_res_list_Mut4.rds")
#homotypic.prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.025*nrow(seurat_obj@meta.data))  ## Assuming 2.5% doublet formation rate - tailor for your dataset
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.025, pK = 0.21, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
#DimPlot(seurat_obj,group.by = c("DF.classifications_0.025_0.06_229","DF.classifications_0.025_0.21_229")) 
WT1_doublet<-seurat_obj

seurat_obj <-WT2 
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
#seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
#seurat_obj <- FindClusters(seurat_obj, resolution = 0.01)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
#optimization of pK. optimized pK value is applied to doubleFinder()
#sweep.res.list_seu <- paramSweep(seurat_obj, PCs = 1:30, sct = FALSE)
#sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
#bcmvn_seu <- find.pK(sweep.stats_seu)
#saveRDS(sweep.res.list_seu,file = "C:/Users/admin/Desktop/sweep_res_list_Mut4.rds")
#homotypic.prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.025*nrow(seurat_obj@meta.data))  ## Assuming 2.5% doublet formation rate - tailor for your dataset
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.025, pK = 0.005, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
#DimPlot(seurat_obj,group.by = c("DF.classifications_0.025_0.06_229","DF.classifications_0.025_0.21_229")) 
WT2_doublet<-seurat_obj

seurat_obj <-Mut3 
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
#seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
#seurat_obj <- FindClusters(seurat_obj, resolution = 0.01)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
#optimization of pK. optimized pK value is applied to doubleFinder()
#sweep.res.list_seu <- paramSweep(seurat_obj, PCs = 1:30, sct = FALSE)
#sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
#bcmvn_seu <- find.pK(sweep.stats_seu)
#saveRDS(sweep.res.list_seu,file = "C:/Users/admin/Desktop/sweep_res_list_Mut4.rds")
#homotypic.prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.025*nrow(seurat_obj@meta.data))  ## Assuming 2.5% doublet formation rate - tailor for your dataset
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.025, pK = 0.06, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
#DimPlot(seurat_obj,group.by = c("DF.classifications_0.025_0.06_229","DF.classifications_0.025_0.21_229")) 
Mut3_doublet<-seurat_obj

seurat_obj <-Mut4 
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
#seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
#seurat_obj <- FindClusters(seurat_obj, resolution = 0.01)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
#optimization of pK. optimized pK value is applied to doubleFinder()
#sweep.res.list_seu <- paramSweep(seurat_obj, PCs = 1:30, sct = FALSE)
#sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
#bcmvn_seu <- find.pK(sweep.stats_seu)
#saveRDS(sweep.res.list_seu,file = "C:/Users/admin/Desktop/sweep_res_list_Mut4.rds")
#homotypic.prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.025*nrow(seurat_obj@meta.data))  ## Assuming 2.5% doublet formation rate - tailor for your dataset
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.025, pK = 0.06, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
#DimPlot(seurat_obj,group.by = c("DF.classifications_0.025_0.06_229","DF.classifications_0.025_0.21_229")) 
Mut4_doublet<-seurat_obj

seurat_obj <-merge(x = WT1_doublet, y = c(WT2_doublet, Mut3_doublet, Mut4_doublet)) #or seperate libraries for appendix
#saveRDS(seurat_obj,file = "C:/Users/admin/Desktop/fulldata_doublet.rds")
doubletmethod2<-seurat_obj
doubletmethod2@meta.data$doublet_merge<-coalesce(doubletmethod2@meta.data$DF.classifications_0.025_0.21_209,
                                                 doubletmethod2@meta.data$DF.classifications_0.025_0.005_325,
                                                 doubletmethod2@meta.data$DF.classifications_0.025_0.06_229,
                                                 doubletmethod2@meta.data$DF.classifications_0.025_0.06_217,)
idx<-match(rownames(fulldata@meta.data),rownames(doubletmethod2@meta.data))
fulldata@meta.data$DF_M2<-doubletmethod2@meta.data[idx,"doublet_merge"]


#####################method2: doublet detection following MAD filter
#MAD filtering
seurat_obj <- WT1 #run following lines 4 times for each library!!!
seurat_obj <- as.SingleCellExperiment(seurat_obj)
mt_genes<-rownames(seurat_obj)[grepl("^mt-", rownames(seurat_obj))]
seurat_obj<-addPerCellQC(seurat_obj,subsets = list(mito = mt_genes))
lib_out <- isOutlier(seurat_obj$sum, nmads = 3, type = "both")
feat_out <- isOutlier(seurat_obj$detected, nmads = 3, type = "both")
discard <- lib_out|feat_out
seurat_obj <- seurat_obj[, !discard]
seurat_obj <- as.Seurat(seurat_obj)
seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj,pattern = "^mt-")
seurat_obj<-subset(seurat_obj, subset=percent.mt <5)
#running reduction
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
#doublet finding
#optimization of pK. optimized pK value is applied to doubleFinder()
#sweep.res.list_seu <- paramSweep(seurat_obj, PCs = 1:30, sct = FALSE)
#sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
#bcmvn_seu <- find.pK(sweep.stats_seu)
#saveRDS(sweep.res.list_seu,file = "C:/Users/admin/Desktop/sweep_res_list_Mut4.rds")
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.025, 
                            pK = 0.19, 
                            nExp = round(0.025*nrow(seurat_obj@meta.data)), 
                            reuse.pANN = NULL, sct = FALSE)
WT1_doublet<-seurat_obj

#MAD filtering
seurat_obj <- WT2 #run following lines 4 times for each library!!!
seurat_obj <- as.SingleCellExperiment(seurat_obj)
mt_genes<-rownames(seurat_obj)[grepl("^mt-", rownames(seurat_obj))]
seurat_obj<-addPerCellQC(seurat_obj,subsets = list(mito = mt_genes))
lib_out <- isOutlier(seurat_obj$sum, nmads = 3, type = "both")
feat_out <- isOutlier(seurat_obj$detected, nmads = 3, type = "both")
discard <- lib_out|feat_out
seurat_obj <- seurat_obj[, !discard]
seurat_obj <- as.Seurat(seurat_obj)
seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj,pattern = "^mt-")
seurat_obj<-subset(seurat_obj, subset=percent.mt <5)
#running reduction
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
#doublet finding
#optimization of pK. optimized pK value is applied to doubleFinder()
#sweep.res.list_seu <- paramSweep(seurat_obj, PCs = 1:30, sct = FALSE)
#sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
#bcmvn_seu <- find.pK(sweep.stats_seu)
#saveRDS(sweep.res.list_seu,file = "C:/Users/admin/Desktop/sweep_res_list_Mut4.rds")
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.025, 
                            pK = 0.01, 
                            nExp = round(0.025*nrow(seurat_obj@meta.data)), 
                            reuse.pANN = NULL, sct = FALSE)
WT2_doublet<-seurat_obj

seurat_obj <- Mut3 #run following lines 4 times for each library!!!
seurat_obj <- as.SingleCellExperiment(seurat_obj)
mt_genes<-rownames(seurat_obj)[grepl("^mt-", rownames(seurat_obj))]
seurat_obj<-addPerCellQC(seurat_obj,subsets = list(mito = mt_genes))
lib_out <- isOutlier(seurat_obj$sum, nmads = 3, type = "both")
feat_out <- isOutlier(seurat_obj$detected, nmads = 3, type = "both")
discard <- lib_out|feat_out
seurat_obj <- seurat_obj[, !discard]
seurat_obj <- as.Seurat(seurat_obj)
seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj,pattern = "^mt-")
seurat_obj<-subset(seurat_obj, subset=percent.mt <5)
#running reduction
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
#doublet finding
#optimization of pK. optimized pK value is applied to doubleFinder()
#sweep.res.list_seu <- paramSweep(seurat_obj, PCs = 1:30, sct = FALSE)
#sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
#bcmvn_seu <- find.pK(sweep.stats_seu)
#saveRDS(sweep.res.list_seu,file = "C:/Users/admin/Desktop/sweep_res_list_Mut4.rds")
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.025, 
                            pK = 0.25, 
                            nExp = round(0.025*nrow(seurat_obj@meta.data)), 
                            reuse.pANN = NULL, sct = FALSE)
Mut3_doublet<-seurat_obj

seurat_obj <- Mut4 #run following lines 4 times for each library!!!
seurat_obj <- as.SingleCellExperiment(seurat_obj)
mt_genes<-rownames(seurat_obj)[grepl("^mt-", rownames(seurat_obj))]
seurat_obj<-addPerCellQC(seurat_obj,subsets = list(mito = mt_genes))
lib_out <- isOutlier(seurat_obj$sum, nmads = 3, type = "both")
feat_out <- isOutlier(seurat_obj$detected, nmads = 3, type = "both")
discard <- lib_out|feat_out
seurat_obj <- seurat_obj[, !discard]
seurat_obj <- as.Seurat(seurat_obj)
seurat_obj[["percent.mt"]]<-PercentageFeatureSet(seurat_obj,pattern = "^mt-")
seurat_obj<-subset(seurat_obj, subset=percent.mt <5)
#running reduction
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
#doublet finding
#optimization of pK. optimized pK value is applied to doubleFinder()
#sweep.res.list_seu <- paramSweep(seurat_obj, PCs = 1:30, sct = FALSE)
#sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
#bcmvn_seu <- find.pK(sweep.stats_seu)
#saveRDS(sweep.res.list_seu,file = "C:/Users/admin/Desktop/sweep_res_list_Mut4.rds")
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.025, 
                            pK = 0.005, 
                            nExp = round(0.025*nrow(seurat_obj@meta.data)), 
                            reuse.pANN = NULL, sct = FALSE)
Mut4_doublet<-seurat_obj

seurat_obj <-merge(x = WT1_doublet, y = c(WT2_doublet, Mut3_doublet, Mut4_doublet)) #or seperate libraries for appendix
doubletmethod1<-seurat_obj
doubletmethod1@meta.data$doublet_merge<-coalesce(doubletmethod1@meta.data$DF.classifications_0.025_0.19_198,
                                                 doubletmethod1@meta.data$DF.classifications_0.025_0.01_308,
                                                 doubletmethod1@meta.data$DF.classifications_0.025_0.25_213,
                                                 doubletmethod1@meta.data$DF.classifications_0.025_0.005_195)
idx<-match(rownames(fulldata@meta.data),rownames(doubletmethod1@meta.data))
fulldata@meta.data$DF_M1<-doubletmethod1@meta.data[idx,"doublet_merge"]


#===================================clustering to detect cluster of habenular cells================
##################method 1: clustering on unfiltered merged data --- Louvain algorithm
seurat_obj <-fulldata
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30,k.param = 20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.009)

idx<-match(rownames(fulldata@meta.data),rownames(seurat_obj@meta.data))
fulldata@meta.data$unfilt_mergeclust<-seurat_obj@meta.data[idx,"seurat_clusters"]


#####################method 2: clustering on unfiltered merged data --- K-means
set.seed(123)
seurat_obj <-fulldata
seurat_obj <- RunPCA(seurat_obj, 
                     features = VariableFeatures(seurat_obj),
                     npcs = 20)
km <- kmeans(seurat_obj@reductions[["pca"]]@cell.embeddings[,1:20],centers=4)
seurat_obj$Kmeans <- km$cluster

idx<-match(rownames(fulldata@meta.data),rownames(seurat_obj@meta.data))
fulldata@meta.data$Kmeans<-seurat_obj@meta.data[idx,"Kmeans"]


########################method 3: clustering on filtered seperate data
seurat_obj <- subset(fulldata, subset = genotype =="WT")
seurat_obj <- subset(seurat_obj, subset = !is.na(miQC_mergeclust))
seurat_obj <- subset(seurat_obj, subset = DF_M1=="Singlet")

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30,k.param = 20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.002)

idx<-match(rownames(fulldata@meta.data),rownames(seurat_obj@meta.data))
fulldata@meta.data$filt_sepclust<-seurat_obj@meta.data[idx,"seurat_clusters"]

seurat_obj <- subset(fulldata, subset = genotype =="mutant")
seurat_obj <- subset(seurat_obj, subset = !is.na(miQC_mergeclust))
seurat_obj <- subset(seurat_obj, subset = DF_M1=="Singlet")

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30,k.param = 20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.006)

idx<-match(rownames(fulldata@meta.data),rownames(seurat_obj@meta.data))
matched<- which(!is.na(idx))
fulldata@meta.data$filt_sepclust[matched] <- seurat_obj@meta.data[idx[matched], "seurat_clusters"]


######################method 4: clustering on unfiltered seperate data --- Leiden algorithm
seurat_obj <- subset(fulldata, subset = genotype =="WT")
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:25,k.param = 10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.0015, algorithm = 4)

idx<-match(rownames(fulldata@meta.data),rownames(seurat_obj@meta.data))
fulldata@meta.data$unfilt_sepclust<-seurat_obj@meta.data[idx,"seurat_clusters"]

seurat_obj <- subset(fulldata, subset = genotype =="mutant")
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20,k.param = 10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.009, algorithm = 4)

idx<-match(rownames(fulldata@meta.data),rownames(seurat_obj@meta.data))
matched<- which(!is.na(idx))
fulldata@meta.data$unfilt_sepclust <- as.character(fulldata@meta.data$unfilt_sepclust)
fulldata@meta.data$unfilt_sepclust[matched] <- seurat_obj@meta.data[idx[matched], "seurat_clusters"]


#####tidying up and saving results
colnames(fulldata@meta.data)[colnames(fulldata@meta.data)=="MAD_seperateclust"] <- "Scater_MAD"
colnames(fulldata@meta.data)[colnames(fulldata@meta.data)=="miQC_mergeclust"] <- "miQC"
fulldata@meta.data$Scater_MAD <- as.character(fulldata@meta.data$Scater_MAD)
fulldata@meta.data$miQC <- as.character(fulldata@meta.data$miQC)
fulldata@meta.data$Scater_MAD[is.na(fulldata@meta.data$Scater_MAD)] <- "poor_quality_cell"
fulldata@meta.data$miQC[is.na(fulldata@meta.data$miQC)] <- "poor_quality_cell"
fulldata@meta.data$filt_sepclust <- as.character(fulldata@meta.data$filt_sepclust)
fulldata@meta.data$filt_sepclust[is.na(fulldata@meta.data$filt_sepclust)] <- "poor_quality_cell"
fulldata@meta.data$DF_M1 <- as.character(fulldata@meta.data$DF_M1)
fulldata@meta.data$DF_M1[is.na(fulldata@meta.data$DF_M1)] <- "poor_quality_cell"
fulldata@meta.data
#saveRDS(fulldata,file="C:/Users/admin/Desktop/preprocessed_fulldata_v3.rds")




#=============================plotting===================================================
fulldata<-readRDS(file="C:/Users/admin/Desktop/preprocessed_fulldata_v3.rds")


#0. a quick look of filtering and clustering result
DimPlot(fulldata,group.by = c("DF_M1","DF_M2","unfilt_mergeclust","Kmeans","unfilt_sepclust","filt_sepclust"))


#1. plot unfilt_mergeclust
reduction <- "umap"
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$unfilt_mergeclust
coords$cell <- rownames(coords)
clusters <- unique(coords$cluster)
# Assign default Seurat palette
library(RColorBrewer)
pal <- brewer.pal(n = length(clusters), name = "Dark2")  # or any palette you like
names(pal) <- clusters
pal["0"] <- "#F8766D"
p1<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

#2. plot Kmeans
reduction <- "umap"
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$Kmeans
coords$cluster <- as.factor(coords$cluster)
coords$cell <- rownames(coords)
clusters <- unique(coords$cluster)
# Assign default Seurat palette
library(RColorBrewer)
pal <- brewer.pal(n = length(clusters), name = "Dark2")  # or any palette you like
names(pal) <- clusters
pal["2"] <- "#F8766D"
p2<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

#3. plot unfilt_sepclust
reduction <- "umap"
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$unfilt_sepclust
coords$cluster <- as.factor(coords$cluster)
coords$cell <- rownames(coords)
clusters <- unique(coords$cluster)
# Assign default Seurat palette
library(RColorBrewer)
pal <- brewer.pal(n = length(clusters), name = "Dark2")  # or any palette you like
names(pal) <- clusters
pal["1"] <- "#F8766D"
p3<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

#4. plot filt_sepclust
reduction <- "umap"
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$filt_sepclust
coords$cell <- rownames(coords)
clusters <- unique(coords$cluster)
# Assign default Seurat palette
library(RColorBrewer)
pal <- brewer.pal(n = length(clusters), name = "Dark2")  # or any palette you like
names(pal) <- clusters
# Override poor-quality cluster
pal["poor_quality_cell"] <- "black"
pal["0"] <- "#F8766D"
p4<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

#ggsave("c:/Users/admin/Desktop/p_unfilt_mergeclust.svg",
#       plot = p1,
#       width = 7, height =7)
#ggsave("c:/Users/admin/Desktop/p_Kmeans.svg",
#       plot = p2,
#       width = 7, height =7)
#ggsave("c:/Users/admin/Desktop/p_unfilt_sepclust.svg",
#       plot = p3,
#       width = 7, height =7)
#ggsave("c:/Users/admin/Desktop/p_filt_sepclust.svg",
#       plot = p4,
#       width = 7, height =7)

#5.plot miQC
reduction <- "umap"
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$miQC
coords$cell <- rownames(coords)
coords <-coords[coords$cluster == "poor_quality_cell", ]
clusters <- unique(coords$cluster)
pal <- c(  "poor_quality_cell" = "red")
p_miQC<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )

#6.plot Scater_MAD
reduction <- "umap"
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$Scater_MAD
coords$cell <- rownames(coords)
coords <-coords[coords$cluster == "poor_quality_cell", ]
clusters <- unique(coords$cluster)
pal <- c("poor_quality_cell" = "blue")
p_MAD<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )

#ggsave("c:/Users/admin/Desktop/p_miQC.svg",
#       plot = p_miQC,
#       width = 7, height =7)
#ggsave("c:/Users/admin/Desktop/p_MAD.svg",
#       plot = p_MAD,
#       width = 7, height =7)



#7.plot DF_M1
reduction <- "umap"
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$DF_M1
coords$cell <- rownames(coords)
clusters <- unique(coords$cluster)
pal <- c("Singlet" = "grey",
         "poor_quality_cell" = "grey",
         "Doublet" = "grey")
p_DFM1_bg<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )
coords <-coords[coords$cluster == "poor_quality_cell", ]
clusters <- unique(coords$cluster)
pal <- c("poor_quality_cell" = "black")
p_DFM1_poor_quality_cell<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$DF_M1
coords$cell <- rownames(coords)
coords <-coords[coords$cluster == "Doublet", ]
clusters <- unique(coords$cluster)
pal <- c("Doublet" = "red")
p_DFM1_Doublet<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )




#8.plot DF_M2
reduction <- "umap"
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$DF_M2
coords$cell <- rownames(coords)
clusters <- unique(coords$cluster)
pal <- c("Singlet" = "grey",
         "Doublet" = "blue")
p_DFM2<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )

#8.plot DF_M2(2 layers image)
reduction <- "umap"
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$DF_M2
coords$cell <- rownames(coords)
clusters <- unique(coords$cluster)
pal <- c("Singlet" = "grey",
         "Doublet" = "gray")
p_DFM2_L1<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )
reduction <- "umap"
coords <- Embeddings(fulldata[[reduction]]) %>% as.data.frame()
coords$cluster <- fulldata$DF_M2
coords$cell <- rownames(coords)
coords <-coords[coords$cluster != "Singlet", ]
clusters <- unique(coords$cluster)
pal <- c("Doublet" = "blue")
p_DFM2_L2_blue<-ggplot(coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )

#ggsave("c:/Users/admin/Desktop/p_DFM1_bg.svg",
#       plot = p_DFM1_bg,
#       width = 7, height =7)
#ggsave("c:/Users/admin/Desktop/p_DFM1_poor_quality_cell.svg",
#       plot = p_DFM1_poor_quality_cell,
#       width = 7, height =7)
#ggsave("c:/Users/admin/Desktop/p_DFM1_Doublet.svg",
#       plot = p_DFM1_Doublet,
#       width = 7, height =7)
#ggsave("c:/Users/admin/Desktop/p_DFM2_L2_blue.svg",
#       plot = p_DFM2_L2_blue,
#       width = 7, height =7)

#9. ploting marker genes
p_Hbmarker<-FeaturePlot(fulldata,c("pou4f1","irx1b","gng8","lmo3"))
p_othermarker<-FeaturePlot(fulldata,c("elavl3","elavl4","epcam","pvalb5"))
#ggsave("c:/Users/admin/Desktop/p_Hbmarker.svg",
#       plot = p_Hbmarker,
#       width = 10, height =10)
#ggsave("c:/Users/admin/Desktop/p_othermarker.svg",
#       plot = p_othermarker,
#       width = 10, height =10)

#==========================stats on habenula-only cluster=======================
fulldata@meta.data$Scater_MAD <- ifelse(fulldata@meta.data$Scater_MAD != "poor_quality_cell", "normal", fulldata@meta.data$Scater_MAD)
fulldata@meta.data$miQC <- ifelse(fulldata@meta.data$miQC != "poor_quality_cell", "normal", fulldata@meta.data$miQC)

Hb<-subset(fulldata,subset = unfilt_mergeclust==0)
table(Hb@meta.data$orig.ident)
library(dplyr)
dplyr::count(Hb@meta.data, orig.ident, DF_M1)
dplyr::count(Hb@meta.data, orig.ident, Scater_MAD)
dplyr::count(Hb@meta.data, orig.ident, miQC)


#==========================plot on merge dataset compare between genotypes/bioreplicates=======================
p_allgenotype<-DimPlot(fulldata,group.by="genotype")
WT<-subset(fulldata,subset=genotype=="WT")
p_WT12<-DimPlot(WT,group.by="orig.ident")
mut<-subset(fulldata,subset=genotype=="mutant")
p_mut12<-DimPlot(mut,group.by="orig.ident")
ggsave("c:/Users/admin/Desktop/p_allgenotype.svg",
       plot = p_allgenotype,
       width = 7, height =7)
ggsave("c:/Users/admin/Desktop/p_WT12.svg",
       plot = p_WT12,
       width = 7, height =7)
ggsave("c:/Users/admin/Desktop/p_mut12.svg",
       plot = p_mut12,
       width = 7, height =7)
#plot the zoomed-in view: plotting olfactory neuron
fulldata <- FindClusters(fulldata, resolution = 0.5)
DimPlot(fulldata,label = TRUE)
of<-subset(fulldata,idents=9)
seurat_obj<-of
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:25,k.param = 8)
seurat_obj <- FindClusters(seurat_obj, resolution = 2)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:25)
p_zoom_of<-DimPlot(seurat_obj,group.by = c("seurat_clusters","genotype"))
#plot the zoomed-in view: plotting olfactory neuron
vHb<-subset(fulldata,idents=2)
seurat_obj<-vHb
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:25,k.param = 8)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:25)
p_zoom_vHb<-DimPlot(seurat_obj,group.by = c("seurat_clusters","orig.ident"))
#plot zoomed in only in WT
vHb_WT<-subset(vHb,subset=genotype=="WT")
seurat_obj<-vHb_WT
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:25,k.param = 8)
seurat_obj <- FindClusters(seurat_obj, resolution = 1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:25)
p_zoom_vHb_WT<-DimPlot(seurat_obj,group.by = c("seurat_clusters","orig.ident"))
#plot zoomed in only in mutant
vHb_mut<-subset(vHb,subset=genotype=="mutant")
seurat_obj<-vHb_mut
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:25,k.param = 8)
seurat_obj <- FindClusters(seurat_obj, resolution = 1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:25)
p_zoom_vHb_mut<-DimPlot(seurat_obj,group.by = c("seurat_clusters","orig.ident"))


