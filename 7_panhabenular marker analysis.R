#- Project 'E:/seurat/seurat_basics' loaded. [renv 1.1.4]
#[Workspace loaded from E:/seurat/seurat_basics/.RData]

library(Seurat)
library(Matrix)  # for sparse matrix support, though Seurat will handle it
library(presto)
#install.packages("svglite")
library(svglite)

#==================================prepare 6dpf telencephalon data========================
# Paths to the files (adjust as needed)
barcode_path1 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523903_6dpf_b1_barcodes.tsv"
feature_path1 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523903_6dpf_b1_features.tsv"
matrix_path1 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523903_6dpf_b1_matrix.mtx"
barcode_path2 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523904_6dpf_b2_barcodes.tsv"
feature_path2 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523904_6dpf_b2_features.tsv"
matrix_path2 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523904_6dpf_b2_matrix.mtx"
barcode_path3 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523905_6dpf_b3_barcodes.tsv"
feature_path3 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523905_6dpf_b3_features.tsv"
matrix_path3 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523905_6dpf_b3_matrix.mtx"
barcode_path4 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523906_6dpf_b4_barcodes.tsv"
feature_path4 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523906_6dpf_b4_features.tsv"
matrix_path4 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523906_6dpf_b4_matrix.mtx"
barcode_path5 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523907_6dpf_b5_barcodes.tsv"
feature_path5 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523907_6dpf_b5_features.tsv"
matrix_path5 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523907_6dpf_b5_matrix.mtx"
barcode_path6 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523908_6dpf_b6_barcodes.tsv"
feature_path6 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523908_6dpf_b6_features.tsv"
matrix_path6 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523908_6dpf_b6_matrix.mtx"
barcode_path7 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523909_6dpf_b7_barcodes.tsv"
feature_path7 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523909_6dpf_b7_features.tsv"
matrix_path7<- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523909_6dpf_b7_matrix.mtx"
barcode_path8 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523910_6dpf_b8_barcodes.tsv"
feature_path8 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523910_6dpf_b8_features.tsv"
matrix_path8 <- "C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523910_6dpf_b8_matrix.mtx"

# Read the data as a sparse matrix
# (ReadMtx is a convenience function in the Seurat package)
expr_mat8 <- ReadMtx(
  mtx = matrix_path8,
  features = feature_path8,
  cells = barcode_path8,
  # optionally, you can set `feature.column = 2` if the features file has gene names in second column
  feature.column = 2
)
# expr_mat should now be a dgCMatrix (genes × cells)
# Create the Seurat object
seurat_obj8 <- CreateSeuratObject(
  counts = expr_mat8,
  project = "6dpf_b8",
  min.cells = 3,
  min.features = 200
)

#saveRDS(seurat_obj1,"C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523903_6dpf_b1.rds")
#saveRDS(seurat_obj2,"C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523904_6dpf_b2.rds")
#saveRDS(seurat_obj3,"C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523905_6dpf_b3.rds")
#saveRDS(seurat_obj4,"C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523906_6dpf_b4.rds")
#saveRDS(seurat_obj5,"C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523907_6dpf_b5.rds")
#saveRDS(seurat_obj6,"C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523908_6dpf_b6.rds")
#saveRDS(seurat_obj7,"C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523909_6dpf_b7.rds")
#saveRDS(seurat_obj8,"C:/Users/admin/Desktop/seurat/6dpf telecephalon/GSM6523910_6dpf_b8.rds")

tel<-merge(x=seurat_obj1,y=c(seurat_obj2,seurat_obj3,seurat_obj4,seurat_obj5,seurat_obj6,seurat_obj7,seurat_obj8))
#saveRDS(tel,"C:/Users/admin/Desktop/seurat/6dpf telecephalon/b1_to_b8_merged.rds")
# (Optional) Add mitochondrial percentage, QC filtering, normalization etc.
#seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
#seurat_obj <- subset(
#  seurat_obj,
#  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10
#)


#=========================processing & find pan-habenula markers=====================
tel<-readRDS("C:/Users/admin/Desktop/seurat/6dpf telecephalon/b1_to_b8_merged.rds")
tel<-JoinLayers(tel)
tel <- NormalizeData(tel)
tel <- FindVariableFeatures(tel, selection.method = "vst", nfeatures = 2000)
tel <- ScaleData(tel)
tel <- RunPCA(tel)
tel <- FindNeighbors(tel, dims = 1:30)
tel <- FindClusters(tel, resolution = 0.005)
tel <- RunUMAP(tel, dims = 1:30)
DimPlot(tel, reduction = "umap", label = TRUE)
FeaturePlot(tel,c("pou4f1","irx1b","gng8","kiss1"))

marker<-FindMarkers(tel, ident.1 = 2)
marker$gene<-rownames(marker)
#write.csv(marker,"C:/Users/admin/Desktop/panhabenular marker.csv")


#=========================plotting===================================================
Hb<-readRDS("C:/Users/admin/Desktop/tidy up code/Hb_allgenotype_withintegration.rds")
WT<-readRDS("C:/Users/admin/Desktop/tidy up code/Hb_WT.rds")
mut<-readRDS("C:/Users/admin/Desktop/tidy up code/Hb_mut.rds")
library(gridExtra)
library(ggplot2)

gene<-rownames(Hb)
marker2<-marker[marker$gene%in%rownames(Hb),]

for (x in 1:100){
  gene<-as.character(marker2[x,"gene"])
  p1<-FeaturePlot(WT,reduction = "umap_cca", gene)
  p2<-FeaturePlot(mut,reduction = "umap_cca", gene)
  p3<-FeaturePlot(tel,gene)
  p<-grid.arrange(p1,p2,p3,ncol=3)
  dir<-file.path("C:/Users/admin/Desktop/panhabenular marker",paste0(gene,".png"))
  ggsave(dir,plot = p,width=21,height=7,dpi=500)
}

p1<-FeaturePlot(WT,reduction = "umap_cca", "CNIH3")
p2<-FeaturePlot(mut,reduction = "umap_cca", "CNIH3")
p3<-FeaturePlot(tel,reduction = "umap", "cnih3")
p<-grid.arrange(p1,p2,p3,ncol=3)


#=============================checking logFC on whole Hb pseudobulk=============================
res_df<-read.csv('C:/Users/admin/Desktop/2025-11-25 correction/DEseq result/WTvsmut/Hb_Deseq_Mut_vs_WT.csv')

pan<-read.csv("C:/Users/admin/Desktop/panhabenular marker.csv")
pan_list<-pan$gene[pan$marker=="marker"]

res_df$panhabenular_marker<-"other"
res_df[res_df$gene %in% pan_list, "panhabenular_marker"]<-"marker"

write.csv(res_df,"C:/Users/admin/Desktop/panhabenular_marker.csv")
