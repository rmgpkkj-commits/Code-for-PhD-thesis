#- Project 'E:/seurat/seurat_basics' loaded. [renv 1.1.4]
#[Workspace loaded from E:/seurat/seurat_basics/.RData]
library(Seurat)
library(presto)
library(svglite)
library(dplyr)
library(ggplot2)
library(svglite)
library(gridExtra)
library(readxl)

#==========================obtaining haebular dataset=================================
Hb<-readRDS("C:/Users/admin/Desktop/tidy up code/Hb_allgenotype_withintegration.rds")
Hb_WT <- subset(Hb, subset = genotype=="WT")
Hb_mut <- subset(Hb, subset = genotype=="mutant")


set.seed(123)
#============================clustering individual WT dataset=========================
seurat_obj <- NormalizeData(Hb_WT)
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst",
                                   assay = "RNA",
                                   nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30,k.param = 20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:35)
Hb_WT<-seurat_obj
Hb_WT@meta.data$cluster<-Hb_WT$seurat_clusters
DimPlot(Hb_WT,label=TRUE, reduction = "umap_cca",group.by = "cluster")
#saveRDS(Hb_WT, file="C:/Users/admin/Desktop/tidy up code/Hb_WT.rds")


#==============================comparing WT dataset with pandey dataset=================
pandey_marker<-read_xlsx("C:/Users/admin/Desktop/tidy up code/pandey_marker.xlsx")
pandey_marker<-pandey_marker[!is.na(pandey_marker$Gene),]
# Define the desired cluster order
cluster_order <- c(0, 1, 10, 5, 9, 6, 7, 4, 2, 11, 12, 8, 3)
# Convert identities to a factor with that specific order
Idents(Hb_WT) <- factor(Idents(Hb_WT), levels = cluster_order)
p_heatmap<-DoHeatmap(Hb_WT,pandey_marker$Gene,assay = "RNA")
#警告信息:
#  In DoHeatmap(Hb_WT, pandey_marker$Gene, assay = "RNA") :
#  The following features were omitted as they were not found in the scale.data slot for the RNA assay: prr15la, gstp1.1, s100z, nmu, dlx5a, gad1b, agrp, smox, galn, rpl26, fxyd6l, drd4a, phc2a, bsnb, mtch2, pth2r, gpr78, c1ql4b, dcdc2b, cyth4a, vav3b, ttna, impact
#sum(grepl("fxyd1",rownames(Hb_WT@assays[["RNA"]])))
#sum(grepl("kiss1rb",rownames(Hb_WT@assays[["RNA"]]@layers[["scale.data"]])))
ggsave("c:/Users/admin/Desktop/p_heatmap.svg",
       plot = p_heatmap,
       width = 7, height =5)


#==========plotting unintegrated merged dataset (showing identity blurriness========================
P1<-DimPlot(Hb,reduction = "pca",group.by = "genotype")
Hb <-RunUMAP(Hb, dims = 1:35, reduction = "pca",reduction.name = "umap_uninteg")
P2<-DimPlot(Hb,reduction = "umap_uninteg",group.by = "genotype")
ggsave("c:/Users/admin/Desktop/mergPCA.svg",
       plot = P1,
       width = 7, height =7)
ggsave("c:/Users/admin/Desktop/mergUMAP.svg",
       plot = P2,
       width = 7, height =7)


#===================plotting cluster number using different algorithm and parameters================
Hb_WT<-FindClusters(Hb_mut, resolution = 0.005, algorithm = 4)
unique(Hb_WT$seurat_clusters)
# Sample data
x <- c(0.001,0.003,0.005,0.01,0.03,0.08,0.1,0.2,0.3,0.5,0.6) #run standard clustering algorithm 
#FindClusters(Hb_WT or Hb_mut, resolution = numbers in x, algorithm = 1 or 4) using either
#algorithm=1 (for Louvaine cluster) or algorithm=4(for Leiden algorithm) and resolution specified in x, note down
#number of resulting clusters manually, and you will obtain the following results:
Lu_WT <- c(1,2,3,4,5,7,8,9,12,13,15)
Le_WT <- c(1,1,3,4,5,7,8,9,11,13,14)
Lu_mut <- c(1,1,1,2,2,5,5,8,10,12,15)
Le_mut <- c(1,1,1,2,2,4,5,6,8,11,11)
#plotting results: 
df <- data.frame(x, Lu_WT, Le_WT, Lu_mut, Le_mut)
df_long <- pivot_longer(df, cols = Lu_WT:Le_mut, names_to = "series", values_to = "y")
line_types <- c("Lu_WT" = "solid", "Le_WT" = "dashed", "Lu_mut" = "solid", "Le_mut" = "dashed")
my_colors <- c("Lu_WT" = "red", "Le_WT" = "red", "Lu_mut" = "blue", "Le_mut" = "blue")
p_clusternumber <- ggplot(df_long, aes(x = x, y = y, color = series, linetype = series)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +           # show points
  scale_linetype_manual(values = line_types) +
  scale_color_manual(values = my_colors) +
  labs(title = "Line Plot with 4 Lines",
       x = "X-axis",
       y = "Y-axis") +
  theme_minimal()
ggsave("c:/Users/admin/Desktop/clusternumber.svg",
       plot = p_clusternumber,
       width = 10, height =5)


#================================Silhouette distance===================================
library(cluster)
#preparing WT and mutant  dataset for silhouette distance measurement
set.seed(123)  # for reproducibility
seurat_obj <- subset(
  Hb_WT,
  cells = sample(colnames(Hb_WT), 12403) #to subsample WT dataset so that mutant and WT cell number are the same
)
seurat_obj <-CreateSeuratObject(seurat_obj[["RNA"]],meta.data = seurat_obj@meta.data)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst",
                                   assay = "RNA",
                                   nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30,k.param = 20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:35)
Sil_WT<-seurat_obj

seurat_obj <-CreateSeuratObject(Hb_mut[["RNA"]],meta.data = Hb_mut@meta.data)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst",
                                   assay = "RNA",
                                   nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30,k.param = 20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:35)
Sil_mut<-seurat_obj

#silhouette distance measurement for WT
Sil_WT <- FindClusters(Sil_WT, resolution = 0.8,algorithm = 4) #try different resolution and algorithm
emb<-Embeddings(Sil_WT,"pca")
clusters<-Sil_WT$seurat_clusters
clusters_numeric<-as.numeric(as.factor(clusters))
dist_matrix<-dist(emb)
sil<-silhouette(clusters_numeric,dist_matrix)
Sil_WT$sil_Leiden_0.8<-sil[,"sil_width"]   #change column accordingly to save results in seurat object
summary(sil)
p3<-ggplot(Sil_WT@meta.data, aes(x = genotype, y = sil_Leiden_0.8, fill = genotype)) + ###change y accordingly
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_classic() +
  labs(y = "Silhouette width", x = "Condition")

#silhouette distance measurement for mutant
Sil_mut <- FindClusters(Sil_mut, resolution = 0.8, algorithm = 4) #try different resolution and algorithm
emb<-Embeddings(Sil_mut,"pca")
clusters<-Sil_mut$seurat_clusters
clusters_numeric<-as.numeric(as.factor(clusters))
dist_matrix<-dist(emb)
sil<-silhouette(clusters_numeric,dist_matrix)
Sil_mut$sil_Leiden_0.8<-sil[,"sil_width"]   #change column accordingly to save results in seurat object
summary(sil)
p4<-ggplot(Sil_mut@meta.data, aes(x = genotype, y = sil_Leiden_0.8, fill = genotype)) + #change y accordingly
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_classic() +
  labs(y = "Silhouette width", x = "Condition")

#saveRDS(Sil_WT,file="C:/Users/admin/Desktop/Sil_WT.rds")
#saveRDS(Sil_mut,file="C:/Users/admin/Desktop/Sil_mut.rds")

#plotting silhouette distance results
library(ggridges)
library(ggplot2)
#Sil_WT<-readRDS(file="C:/Users/admin/Desktop/Sil_WT.rds")
#Sil_mut<-readRDS(file="C:/Users/admin/Desktop/Sil_mut.rds")

df_WT<-Sil_WT@meta.data[,c("sil_Louvaine_0.05","sil_Louvaine_0.1","sil_Louvaine_0.3","sil_Louvaine_0.8","sil_Louvaine_1.5",
                           "sil_Leiden_0.05","sil_Leiden_0.1","sil_Leiden_0.3","sil_Leiden_0.8","sil_Leiden_1.5","genotype")]
df_mut<-Sil_mut@meta.data[,c("sil_Louvaine_0.05","sil_Louvaine_0.1","sil_Louvaine_0.3","sil_Louvaine_0.8","sil_Louvaine_1.5",
                             "sil_Leiden_0.05","sil_Leiden_0.1","sil_Leiden_0.3","sil_Leiden_0.8","sil_Leiden_1.5","genotype")]
df<-rbind(df_WT,df_mut)
library(tidyverse)
df <- df %>%
  pivot_longer(
    cols = starts_with("sil_"),              # all silhouette columns
    names_to = "metric",                     # new column: which sil_ column
    values_to = "silhouette_width"           # new column: the values
  )

ggplot(df, aes(x = silhouette_width, y = metric, fill = genotype)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
  theme_classic() +
  labs(x = "Silhouette width", y = "Resolution / Method", fill = "Genotype")

means_df <- df %>%
  group_by(metric, genotype) %>%
  summarize(mean_sil = mean(silhouette_width, na.rm = TRUE))
P_ridge<-ggplot(df, aes(x = silhouette_width, y = metric, fill = genotype)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
  geom_segment(
    data = means_df,
    aes(
      x = mean_sil,
      xend = mean_sil,
      y = as.numeric(factor(metric)) - 0.2,  # small offset down
      yend = as.numeric(factor(metric)) + 0.2, # small offset up
      color = genotype
    ),
    linewidth = 1.1
  ) +
  scale_fill_manual(values = c("WT" = "#1f78b4", "Mut" = "#e31a1c")) +
  scale_color_manual(values = c("WT" = "black", "Mut" = "red")) +
  theme_classic() +
  labs(x = "Silhouette width", y = "Resolution / Method")

ggsave("c:/Users/admin/Desktop/silhouette.svg",
       plot = P_ridge,
       width = 10, height =5)




#Hb_WT<-readRDS( file="C:/Users/admin/Desktop/Hb_WT.rds")
Hb_WT_integ <- subset(Hb, subset = genotype=="WT")
Hb_mut_integ <- subset(Hb, subset = genotype=="mutant")
Hb_WT[["umap_integ"]]<-Hb_WT_integ[["umap_cca"]]
Hb_mut[["umap_integ"]]<-Hb_mut_integ[["umap_cca"]]
#========================mutant clustering (independent/label transfer)=====================
#Hb_mut <- subset(Hb, subset = genotype=="mutant")
#annotating mutant library using label transfer
seurat_obj <- FindVariableFeatures(Hb_mut)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:35)
anchor <-FindTransferAnchors(reference = Hb_WT,
                             query = seurat_obj,
                             dims = 1:30,
                             reference.reduction = "pca")
predictions <- TransferData(anchorset = anchor, 
                            refdata = Hb_WT$seurat_clusters,
                            dims = 1:30)
seurat_obj<- AddMetaData(seurat_obj,metadata = predictions)
Hb_mut<-seurat_obj
Hb_mut$predicted.id <-factor(Hb_mut$predicted.id, levels = as.character(0:12))
DimPlot(Hb_mut,reduction = "umap_cca",group.by = "predicted.id", label = TRUE)

#clustering mutant library using standard work flow
Hb_mut <- FindNeighbors(Hb_mut, dims = 1:30,k.param = 20)
Hb_mut <- FindClusters(Hb_mut, resolution = 0.4,algorithm = 1)
Hb_mut$Louvaine_0.4<-Hb_mut$seurat_clusters
Hb_mut <-FindSubCluster(Hb_mut, cluster = "6", graph.name= "RNA_snn",resolution = 0.08)
Hb_mut$seurat_clusters<-Hb_mut$sub.cluster
Idents(Hb_mut)<-"seurat_clusters"
Hb_mut <-FindSubCluster(Hb_mut, cluster = "4", graph.name= "RNA_snn",resolution = 0.08)
Hb_mut@meta.data$cluster<-Hb_mut$sub.cluster
DimPlot(Hb_mut, reduction = "umap_cca",group.by="cluster",label=TRUE)

#saveRDS(Hb_mut, file="C:/Users/admin/Desktop/tidy up code/Hb_mut.rds")

#plotting clustering and label transfer results
#unique(ggplot_build(p)$data[[1]]$colour)
cluster_colors <- c(
  "0" = "#8CAB00",
  "1" = "#F8766D",
  "2" = "#E18A00",
  "3" = "#00ACFC",
  "4_0" = "#24B700",
  "4_1"="#00BBDA",
  "5" = "#8B93FF",
  "6_0" = "#BE9C00",
  "6_1" = "#F962DD",
  "7" = "#F781BF",
  "8" = "#00BE70",
  "9" = "#00C1AB"
)
p_clust<-DimPlot(Hb_WT, reduction="umap_cca",group.by = "cluster",label =FALSE)+
  DimPlot(Hb_mut, reduction = "umap_cca",group.by = "predicted.id",label=FALSE)+
  DimPlot(Hb_mut, reduction = "umap_cca",group.by = "cluster",label=FALSE, cols = cluster_colors)

#plotting prediction score
p_predictionscore <- FeaturePlot(Hb_mut, reduction= "umap_cca", 
                                 features = "prediction.score.max", 
                                 cols = c("blue", "orange"))

ggsave("c:/Users/admin/Desktop/cluster.svg",
       plot = p_clust,
       width = 24, height =8)
ggsave("c:/Users/admin/Desktop/predictionscore.svg",
       plot = p_predictionscore,
       width = 8, height =8)

#ploting Hb06 and vHb marker genes (to show some clusters can still be aligned across genotypes)
p_wnt11_prkcq_gng8_wt <- FeaturePlot(Hb_WT,c("wnt11","prkcq","gng8"),reduction="umap_integ")
p_wnt11_prkcq_gng8_mut <- FeaturePlot(Hb_mut,c("wnt11","prkcq","gng8"),reduction="umap_integ")
p_pou3f1_pnoca<-FeaturePlot(Hb_WT,c("pou3f1","ENSDARG00000025024"),blend = TRUE,reduction="umap_integ")[[3]]+
  FeaturePlot(Hb_mut,c("pou3f1","ENSDARG00000025024"),blend = TRUE,reduction="umap_integ")[[3]]
p_foxa1_gng8<-FeaturePlot(Hb_WT,c("foxa1","gng8"),blend = TRUE,reduction="umap_integ")[[3]]+
  FeaturePlot(Hb_mut,c("foxa1","gng8"),blend = TRUE,reduction="umap_integ")[[3]]
p_color<-FeaturePlot(Hb_WT,c("pou3f1","ENSDARG00000025024"),blend = TRUE,reduction="umap_integ")[[4]]
ggsave("c:/Users/admin/Desktop/wnt11_prkcq_gng8_wt.svg",
       plot = p_wnt11_prkcq_gng8_wt,
       width = 16, height =16)
ggsave("c:/Users/admin/Desktop/wnt11_prkcq_gng8_mut.svg",
       plot = p_wnt11_prkcq_gng8_mut,
       width = 16, height =16)
ggsave("c:/Users/admin/Desktop/pou3f1_pnoca.svg",
       plot = p_pou3f1_pnoca,
       width = 16, height =8)
ggsave("c:/Users/admin/Desktop/foxa1_gng8.svg",
       plot = p_foxa1_gng8,
       width = 16, height =8)
ggsave("c:/Users/admin/Desktop/colorthreshold.svg",
       plot = p_color,
       width = 5, height =5)



#========================checking clustering stability================================
Hb_WT<-readRDS(file="C:/Users/admin/Desktop/Hb_WT.rds")
Hb_mut<-readRDS(file="C:/Users/admin/Desktop/Hb_mut.rds")
WT_0.0625<-subset(Hb_WT,cells = sample(colnames(Hb_WT), round(0.0625*16779)))
WT_0.125<-subset(Hb_WT,cells = sample(colnames(Hb_WT), round(0.125*16779)))
WT_0.25<-subset(Hb_WT,cells = sample(colnames(Hb_WT), round(0.25*16779)))
WT_0.375<-subset(Hb_WT,cells = sample(colnames(Hb_WT), round(0.375*16779)))
WT_0.5<-subset(Hb_WT,cells = sample(colnames(Hb_WT), round(0.5*16779)))
WT_0.625<-subset(Hb_WT,cells = sample(colnames(Hb_WT), round(0.625*16779)))
WT_0.75<-subset(Hb_WT,cells = sample(colnames(Hb_WT), round(0.75*16779)))
WT_0.875<-subset(Hb_WT,cells = sample(colnames(Hb_WT), round(0.875*16779)))
mut_0.0625<-subset(Hb_mut,cells = sample(colnames(Hb_mut), round(0.0625*12403)))
mut_0.125<-subset(Hb_mut,cells = sample(colnames(Hb_mut), round(0.125*12403)))
mut_0.25<-subset(Hb_mut,cells = sample(colnames(Hb_mut), round(0.25*12403)))
mut_0.375<-subset(Hb_mut,cells = sample(colnames(Hb_mut), round(0.375*12403)))
mut_0.5<-subset(Hb_mut,cells = sample(colnames(Hb_mut), round(0.5*12403)))
mut_0.625<-subset(Hb_mut,cells = sample(colnames(Hb_mut), round(0.625*12403)))
mut_0.75<-subset(Hb_mut,cells = sample(colnames(Hb_mut), round(0.75*12403)))
mut_0.875<-subset(Hb_mut,cells = sample(colnames(Hb_mut), round(0.875*12403)))

#run following clustering for each subsampled dataset above and mannually note down cluster number
seurat_obj <-WT_0.0625   #change WT_0.625 to each of the above subsampled dataset
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst",
                                   assay = "RNA",
                                   nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30,k.param = 20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
unique(seurat_obj$seurat_clusters)

# plotting subsampling results
x <- c(0.0625,0.125,0.25,0.5,1) # % of dataset subsampling
n_WT <- c(10,11,13,13,13) # resulting cluster number noted down from above process for WT
n_mut <- c(5,8,10,10,10) # resulting cluster number noted down from above process for mutant
df <- data.frame(x, n_WT, n_mut)
df_long <- pivot_longer(df, cols = n_WT:n_mut, names_to = "series", values_to = "y")
my_colors <- c("n_WT" = "red", "n_mut" = "blue")
p_subset_n <- ggplot(df_long, aes(x = x, y = y, color = series, linetype = series)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +           # show points
  scale_color_manual(values = my_colors) +
  labs(title = "Line Plot with 4 Lines",
       x = "X-axis",
       y = "Y-axis") +
  theme_minimal()
ggsave("c:/Users/admin/Desktop/subset_n.svg",
       plot = p_subset_n,
       width = 10, height =5)


#=========================hierachical clustering==================================
WT<-Hb_WT

WT<-FindClusters(WT,resolution = 0.002, graph.name = "RNA_snn")
WT$L1<-WT$seurat_clusters
WT<-FindClusters(WT,resolution = 0.003, graph.name = "RNA_snn")
WT$L2<-WT$seurat_clusters
WT<-FindClusters(WT,resolution = 0.5, graph.name = "RNA_snn")
WT$clusters<-WT$seurat_clusters
DimPlot(WT, reduction = "umap_cca", label = TRUE)

p1<-DimPlot(WT, reduction = "umap_cca", label = FALSE, group.by = "L1")
p2<-DimPlot(WT, reduction = "umap_cca", label = FALSE, group.by = "L2")
p3<-DimPlot(WT, reduction = "umap_cca", label = FALSE, group.by = "clusters")
p4<-FeaturePlot(WT,reduction = "umap_cca", c("gng8","kiss1","kctd12.1","kctd12.2"))

ggsave("C:/Users/admin/Desktop/L1cluster.png",plot = p1,width=8,height=8,dpi=700)
ggsave("C:/Users/admin/Desktop/L2cluster.png",plot = p2,width=8,height=8,dpi=700)
ggsave("C:/Users/admin/Desktop/finalcluster.png",plot = p3,width=8,height=8,dpi=700)
ggsave("C:/Users/admin/Desktop/canonical marker.png",plot = p4,width=16,height=16,dpi=700)



