#- Project 'E:/seurat/seurat_harmony' loaded. [renv 1.1.4]
#install.packages('harmony')
library(harmony)
library(Seurat)
library(presto)
library(svglite)
library(dplyr)
library(ggplot2)
library(svglite)
library(gridExtra)
library(readxl)

#==========================obtaining haebular dataset=================================
fulldata<-readRDS(file="C:/Users/admin/Desktop/tidy up code/preprocessed_fulldata_v3.rds")
Hb<-subset(fulldata,subset = unfilt_mergeclust==0)

Hb_WT <- subset(Hb, subset = genotype=="WT")
Hb_mut <- subset(Hb, subset = genotype=="mutant")


#====================================integration====================================
#####################1. cca integration##################
Hb[["RNA"]] <- split(Hb[["RNA"]],f=Hb$genotype)
Hb <- NormalizeData(Hb)
Hb <- FindVariableFeatures(Hb)
Hb <- ScaleData(Hb)
Hb <- RunPCA(Hb)
Hb <- IntegrateLayers(object = Hb,
                      method = CCAIntegration,
                      orig.reduction = "pca",
                      new.reduction = "integrated.cca")
Hb[["RNA"]] <- JoinLayers(Hb[["RNA"]])
Hb <- RunUMAP(Hb, dims = 1:35, 
              reduction = "integrated.cca",
              reduction.name = "umap_cca")

####################2. harmony integration##################
Hb<-RunHarmony(Hb, "genotype")
Hb<-RunUMAP(Hb, dims=1:30, reduction = "harmony", reduction.name = "umap_harmony")

####################3. saving results
#saveRDS(Hb,"C:/Users/admin/Desktop/tidy up code/Hb_allgenotype_withintegration.rds")


#========================================running DAseq======================================
#- Project 'E:/seurat/DAseq_seurat' loaded. [renv 1.1.4]
library(Seurat)
library(DAseq)
library(ggplot2)
#Hb<-readRDS("C:/Users/admin/Desktop/tidy up code/Hb_allgenotype_withintegration.rds")

#prearing for pca and umap embeddings
cell.labels<-Hb@meta.data$genotype
labels.1<-"WT"
labels.2<-"mutant"
plot.embedding1<-Hb@reductions[["umap_cca"]]@cell.embeddings
X_cca<-Hb@reductions[["integrated.cca"]]@cell.embeddings
X_harmony<-Hb@reductions[["harmony"]]@cell.embeddings

#running DAseq using integrated_CCA
da_cells_CCA <- getDAcells(
  X = X_cca,
  cell.labels = cell.labels,
  labels.1 = labels.1,
  labels.2 = labels.2,
  k.vector = seq(50,500,50),
  plot.embedding = plot.embedding1
)
#saveRDS(da_cells_CCA, file = "C:/Users/admin/Desktop/tidy up code/result_DAseq_CCA.rds")

#running DAseq using harmony
da_cells_harmony <- getDAcells(
  X = X_harmony,
  cell.labels = cell.labels,
  labels.1 = labels.1,
  labels.2 = labels.2,
  k.vector = seq(50,500,50),
  plot.embedding = plot.embedding1
)
#saveRDS(da_cells_harmony, file = "C:/Users/admin/Desktop/tidy up code/result_DAseq_harmony.rds")


#======================plotting DAseq and integration results=======================================
P1<-DimPlot(Hb,reduction = "umap_cca",group.by = "genotype")
P2<-DimPlot(Hb,reduction = "umap_harmony",group.by = "genotype")

p1<-da_cells_CCA$pred.plot+
  scale_color_gradient2(low="red",mid="yellow",high="blue")
p2<-da_cells_harmony$pred.plot+
  scale_color_gradient2(low="red",mid="yellow",high="blue")

#plot_data<-da_cells_CCA$pred.plot$data
#plot_data[,c("UMAP_1","UMAP_2")]<-plot.embedding1[rownames(plot_data), ]
#ggplot(plot_data,aes(x=UMAP_1,y=UMAP_2,color=Score))+
#  geom_point(size=0.5)+
#  scale_color_gradient2(low="red",mid="yellow",high="blue")

library(svglite)
ggsave("c:/Users/admin/Desktop/cca_DAseq.svg",
       plot = p1,
       width = 8, height =8)
ggsave("c:/Users/admin/Desktop/harmony_DAseq.svg",
       plot = p2,
       width = 8, height =8)
ggsave("c:/Users/admin/Desktop/cca_genotype.svg",
       plot = P1,
       width = 8, height =8)
ggsave("c:/Users/admin/Desktop/harmony_genotype.svg",
       plot = P2,
       width = 8, height =8)


#=========================identify Hb09_0============================================
Hb<-FindNeighbors(Hb,reduction = "integrated.cca", dims = 1:30)
Hb<-FindClusters(Hb,resolution = 1.0, graph.name = "RNA_snn")
#DimPlot(Hb,reduction = "umap_cca",group.by = "seurat_clusters",label = TRUE)
Hb$cca_cluster<-Hb$seurat_clusters
Hb@meta.data[!(Hb$cca_cluster%in%c(6,9)),"cca_cluster"]<-NA
P3<-DimPlot(Hb,reduction = "umap_cca",group.by = "cca_cluster")
#counting WT and mutant cell number in Hb09_0
library(dplyr)
summarize(group_by(Hb@meta.data,genotype,cca_cluster), counts=n())
summarize(group_by(Hb@meta.data,genotype,harmony_cluster), counts=n())

Hb<-FindNeighbors(Hb,reduction = "harmony", dims = 1:30)
Hb<-FindClusters(Hb,resolution = 1.5, graph.name = "RNA_snn",algorithm = 1)
#DimPlot(Hb,reduction = "umap_cca",group.by = "seurat_clusters",label = TRUE)
Hb$harmony_cluster<-Hb$seurat_clusters
Hb@meta.data[!(Hb$harmony_cluster%in%c(15,7,18,23)),"harmony_cluster"]<-NA
P4<-DimPlot(Hb,reduction = "umap_cca",group.by = "harmony_cluster")
#counting WT and mutant cell number in Hb09_0
library(dplyr)
summarize(group_by(Hb@meta.data,genotype,cca_cluster), counts=n())
summarize(group_by(Hb@meta.data,genotype,harmony_cluster), counts=n())

#saving resulting
library(svglite)
ggsave("c:/Users/admin/Desktop/cca_cluster.svg",
       plot = P3,
       width = 8, height =8)
ggsave("c:/Users/admin/Desktop/harmony_cluster.svg",
       plot = P4,
       width = 8, height =8)


