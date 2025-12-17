#- Project 'E:/seurat/DESeq2' loaded. [renv 1.1.4]
library(apeglm)
library(DESeq2)
library(Seurat)
library(presto)
library(org.Dr.eg.db)
library(ggplot2)

#======================grabing list of marker genes======================================
#marker genes for Hb09_0
Hb<-readRDS("C:/Users/admin/Desktop/tidy up code/Hb_allgenotype_withintegration.rds")
Hb<-FindNeighbors(Hb,reduction = "integrated.cca", dims = 1:30)
Hb<-FindClusters(Hb,resolution = 1.0, graph.name = "RNA_snn")
Hb$cca_cluster<-Hb$seurat_clusters
Hb@meta.data[!(Hb$cca_cluster%in%c(6,9)),"cca_cluster"]<-NA

WT<-subset(Hb, subset=genotype=="WT")
Idents(WT)<-"cca_cluster"
Hb09_0<-FindMarkers(WT,ident.1 = "6",assay = "RNA")
Hb09_0$gene<-rownames(Hb09_0)
markerHb09_0<-Hb09_0[Hb09_0$p_val_adj<1e-200,]
markerHb09_0<-markerHb09_0[markerHb09_0$avg_log2FC>0,]
markerHb09_0<-rownames(markerHb09_0)

#marker genes for other clusters #Hb01(0) Hb02(1) Hb05(9) Hb06(6) Hb07(7) Hb09(2) Hb10(11) Hb15(3)
#Hb04(5) Hb08(4) Hb13(8) Hb03(10) Hb11 (12)
WT<-readRDS("C:/Users/admin/Desktop/tidy up code/Hb_WT.rds")
mut<-readRDS("C:/Users/admin/Desktop/tidy up code/Hb_mut.rds")
Idents(WT)<-"seurat_clusters"
Hb01<-FindMarkers(WT,ident.1 = "0",assay = "RNA")
markerHb01<-Hb01[Hb01$p_val_adj<1e-100,]
markerHb01<-markerHb01[markerHb01$avg_log2FC>0,]
markerHb01<-rownames(markerHb01)

Hb02<-FindMarkers(WT,ident.1 = "1",assay = "RNA")
markerHb02<-Hb02[Hb02$p_val_adj<1e-100,]
markerHb02<-markerHb02[markerHb02$avg_log2FC>0,]
markerHb02<-rownames(markerHb02)

Hb09<-FindMarkers(WT,ident.1 = "2",assay = "RNA")
markerHb09<-Hb09[Hb09$p_val_adj<1e-100,]
markerHb09<-markerHb09[markerHb09$avg_log2FC>0,]
markerHb09<-rownames(markerHb09)

Hb15<-FindMarkers(WT,ident.1 = "3",assay = "RNA")
markerHb15<-Hb15[Hb15$p_val_adj<1e-100,]
markerHb15<-markerHb15[markerHb15$avg_log2FC>0,]
markerHb15<-rownames(markerHb15)

Hb08<-FindMarkers(WT,ident.1 = "4",assay = "RNA")
markerHb08<-Hb08[Hb08$p_val_adj<1e-100,]
markerHb08<-markerHb08[markerHb08$avg_log2FC>0,]
markerHb08<-rownames(markerHb08)

Hb04<-FindMarkers(WT,ident.1 = "5",assay = "RNA")
markerHb04<-Hb04[Hb04$p_val_adj<1e-100,]
markerHb04<-markerHb04[markerHb04$avg_log2FC>0,]
markerHb04<-rownames(markerHb04)

Hb06<-FindMarkers(WT,ident.1 = "6",assay = "RNA")
markerHb06<-Hb06[Hb06$p_val_adj<1e-100,]
markerHb06<-markerHb06[markerHb06$avg_log2FC>0,]
markerHb06<-rownames(markerHb06)

Hb07<-FindMarkers(WT,ident.1 = "7",assay = "RNA")
markerHb07<-Hb07[Hb07$p_val_adj<1e-100,]
markerHb07<-markerHb07[markerHb07$avg_log2FC>0,]
markerHb07<-rownames(markerHb07)

Hb13<-FindMarkers(WT,ident.1 = "8",assay = "RNA")
markerHb13<-Hb13[Hb13$p_val_adj<1e-100,]
markerHb13<-markerHb13[markerHb13$avg_log2FC>0,]
markerHb13<-rownames(markerHb13)

Hb05<-FindMarkers(WT,ident.1 = "9",assay = "RNA")
markerHb05<-Hb05[Hb05$p_val_adj<1e-100,]
markerHb05<-markerHb05[markerHb05$avg_log2FC>0,]
markerHb05<-rownames(markerHb05)

Hb03<-FindMarkers(WT,ident.1 = "10",assay = "RNA")
markerHb03<-Hb03[Hb03$p_val_adj<1e-100,]
markerHb03<-markerHb03[markerHb03$avg_log2FC>0,]
markerHb03<-rownames(markerHb03)

Hb10<-FindMarkers(WT,ident.1 = "11",assay = "RNA")
markerHb10<-Hb10[Hb10$p_val_adj<1e-100,]
markerHb10<-markerHb10[markerHb10$avg_log2FC>0,]
markerHb10<-rownames(markerHb10)

Hb11<-FindMarkers(WT,ident.1 = "12",assay = "RNA")
markerHb11<-Hb11[Hb11$p_val_adj<0.0005,]
markerHb11<-markerHb11[markerHb11$avg_log2FC>0,]
markerHb11<-rownames(markerHb11)


#====================fgsea on marker genes using whole Hb psuedobulk============================
#BiocManager::install("fgsea") 
library(fgsea)
res_df<-read.csv('C:/Users/admin/Desktop/2025-11-25 correction/DEseq result/WTvsmut/Hb_Deseq_Mut_vs_WT.csv')

res_df<-res_df[order(res_df$log2FoldChange),]   #order the result by pvalue (optional) or log2FoldChange
gene_list <- res_df$log2FoldChange  # or log2FoldChange
names(gene_list)<-res_df$gene  #names(gene_list)<-res_df$X

pHb15<-plotEnrichment(markerHb15, gene_list)
fgseaRes_Hb15 <- fgsea(
  pathways = list(MyGenes = markerHb15),
  stats = gene_list,
  nperm = 10000
)
fgseaRes_Hb15
pHb15

fgseaRes<-rbind(fgseaRes_Hb01,fgseaRes_Hb02,fgseaRes_Hb03,fgseaRes_Hb04,fgseaRes_Hb05,
                fgseaRes_Hb06,fgseaRes_Hb07,fgseaRes_Hb08,fgseaRes_Hb09,fgseaRes_Hb10,
                fgseaRes_Hb11,fgseaRes_Hb13,fgseaRes_Hb15,fgseaRes_Hb09_0)
fgseaRes$cluster<-c("Hb01","Hb02","Hb03","Hb04","Hb05","Hb06","Hb07","Hb08","Hb09",
                    "Hb10","Hb11","Hb13","Hb15","Hb09_0") 
fgseaRes<-as.data.frame(fgseaRes)
fgseaRes_export <- fgseaRes
fgseaRes_export$leadingEdge <- sapply(
  fgseaRes$leadingEdge,
  function(x) {
    if (length(x) == 0) return(NA)      # empty list
    paste(x, collapse = "/")            # turn into "gene1/gene2/gene3"
  }
)
write.csv(
  fgseaRes_export,
  "C:/Users/admin/Desktop/feseaRes_wholeHb_WTvsmut.csv",
  row.names = FALSE
)

# List of your plots
plots <- list(pHb09_0, pHb01, pHb02, pHb03, pHb04, pHb05, pHb06, pHb07, pHb08, pHb09, pHb10, pHb11, pHb13, pHb15)
plotname<- list("Hb09_0", "Hb01", "Hb02", "Hb03", "Hb04", "Hb05", "Hb06", "Hb07", "Hb08", "Hb09", "Hb10", "Hb11", "Hb13", "Hb15")
# Loop over them and save
for(i in seq_along(plots)) {
  svg_filename <- paste0("C:/Users/admin/Desktop/plot_", plotname[i], ".svg")  # p0 -> plot_0.svg
  svg(svg_filename, width = 8, height = 6)
  print(plots[[i]])   # important to print ggplot/Seurat objects
  dev.off()
}


#======================fgesa on marker genes using cluster-specific pseudobulk================
res_df<-read.csv('C:/Users/admin/Desktop/2025-11-25 correction/DEseq result/WTvsmut/Hb02_Deseq_Mut_vs_WT.csv') #change here to load each cluster-specific pseudobulk dataset
res_df<-res_df[order(res_df$log2FoldChange),]   #order the result by pvalue (optional) or log2FoldChange
gene_list <- res_df$log2FoldChange  # or log2FoldChange
names(gene_list)<-res_df$gene

pHb02<-plotEnrichment(markerHb02, gene_list)
fgseaResHb02 <- fgsea(
  pathways = list(MyGenes = markerHb02),
  stats = gene_list,
  nperm = 10000
)
fgseaResHb02
pHb02

fgseaRes_aligned<-rbind(fgseaResHb02,fgseaResHb04,fgseaResHb05,fgseaResHb06,fgseaResHb07,fgseaResHb08,fgseaResHb10,fgseaResHb15)
fgseaRes_aligned$cluster<-c("Hb02","Hb04","Hb05","Hb06","Hb07","Hb08","Hb10","Hb15") 
fgseaRes_aligned_export <- fgseaRes_aligned
fgseaRes_aligned_export$leadingEdge <- sapply(
  fgseaRes_aligned$leadingEdge,
  function(x) {
    if (length(x) == 0) return(NA)      # empty list
    paste(x, collapse = "/")            # turn into "gene1/gene2/gene3"
  }
)
write.csv(
  fgseaRes_aligned_export,
  "C:/Users/admin/Desktop/aligned_feseaRes_mutvsWT.csv",
  row.names = FALSE
)

plots <- list(pHb02, pHb04, pHb05, pHb06, pHb07, pHb08, pHb10, pHb15)
plotname<- list("Hb02", "Hb04", "Hb05", "Hb06", "Hb07", "Hb08", "Hb10","Hb15")
# Loop over them and save
for(i in seq_along(plots)) {
  svg_filename <- paste0("C:/Users/admin/Desktop/plot_", plotname[i], ".svg")  # p0 -> plot_0.svg
  svg(svg_filename, width = 8, height = 6)
  print(plots[[i]])   # important to print ggplot/Seurat objects
  dev.off()
}


#=========================volcano plot for Hb09_0================================
res_df<-read.csv('C:/Users/admin/Desktop/2025-11-25 correction/DEseq result/WTvsmut/Hb_Deseq_Mut_vs_WT.csv')
res_df$marker<-ifelse(res_df$gene %in% markerHb09_0, "marker", "other")

marker2<-Hb09_0[Hb09_0$avg_log2FC>0,]
marker2<-marker2[marker2$p_val_adj<1e-200,]
marker2 <- marker2[order(marker2$avg_log2FC, decreasing = TRUE), ]
marker2$rank<-c(1:128)
idx<-match(marker2$gene,res_df$gene)
res_df$markerrank<-"other"
res_df$markerrank[idx]<-marker2$rank
res_df2<-res_df[res_df$marker=="marker",]
res_df2$markerrank <- as.numeric(res_df2$markerrank)


# Combine both datasets
combined_y <- c(res_df$padj, res_df2$padj)
combined_y <- combined_y[!is.na(combined_y) & combined_y > 0]
# Compute limits
y_limits <- c(0, max(-log10(combined_y)))
x_limits <- range(c(res_df$log2FoldChange, res_df2$log2FoldChange), na.rm = TRUE)

p1<-ggplot(res_df, aes(x=log2FoldChange,y=-log10(padj),color=marker))+
  geom_point(alpha=0.7, size=1.8)+
  scale_color_manual(values=c("marker"="red", "other"="grey90"))+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+
  geom_hline(yintercept = -log10(0.001), linetype ="dashed")+
  expand_limits(y=0)+
  coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off")+
  theme_minimal(base_size = 14)

p2<-ggplot(res_df2, aes(x = log2FoldChange,
                        y = -log10(padj),
                        color = markerrank)) +
  geom_point(alpha = 1.0, size = 2.3) +
  scale_color_gradientn(
    colours = c("red", "orange", "blue"),values = c(0,0.3,1),
    na.value = "gray90"
  ) +
  geom_vline(xintercept = c(-1,1),linetype="dashed")+
  geom_hline(yintercept = -log10(0.001), linetype ="dashed")+
  #  geom_vline(xintercept = 2.5,linetype="dashed",color="blue")+
  #  geom_hline(yintercept = -log10(1e-100), linetype ="dashed",color="blue")+
  expand_limits(y=0)+
  coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off")+
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot",
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~P~value),
    color = "marker Log2FC"
  ) 


#res_df[is.na(res_df$log2FoldChange) | is.na(res_df$padj),]$marker


library(svglite)
ggsave("c:/Users/admin/Desktop/volcano plot total.svg",
       plot = p1,
       width = 16, height =8)
ggsave("c:/Users/admin/Desktop/volcano plot Hb09_0.svg",
       plot = p2,
       width = 16, height =8)


#finding the 6 specific Hb09_0 marker
res_df2$Hb09_0<-"other"
res_df2$Hb09_0[res_df2$gene=="ENSDARG00000104624"]<-"marker"
ggplot(res_df2, aes(x=log2FoldChange,y=-log10(padj),color=Hb09_0))+
  geom_point(alpha=0.7, size=1.8)+
  scale_color_manual(values=c("marker"="red", "other"="black"))+
  geom_vline(xintercept = c(-1,1),linetype="dashed")+
  geom_hline(yintercept = -log10(0.001), linetype ="dashed")+
  expand_limits(y=0)+
  coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off")+
  theme_minimal(base_size = 14)

#=========================finding genes seperating Hb09_0 and Hb09_1====================
WT<-subset(Hb, subset=genotype=="WT")
mut<-subset(Hb, subset=genotype=="mutant")
Idents(WT)<-"cca_cluster"
marker<-FindMarkers(WT,ident.1 = 6)
marker_filt<-marker[marker$pct.2<0.05,]
marker_filt<-marker_filt[marker_filt$pct.1>0.15,]
marker_filt$gene<-rownames(marker_filt)
marker_filt$absent_in_mutant<-TRUE
marker_filt$absent_in_mutant[
  rownames(marker_filt) %in% 
    c("si:dkey-183i3.6","tac1","chrna4b")] <- FALSE

#marker_filt2<-marker[marker$pct.1>10*marker$pct.2,]
marker_filt3<-marker[marker$pct.1>7*marker$pct.2,]
marker_filt3$gene<-rownames(marker_filt3)
marker_filt3$absent_in_mutant<-TRUE
marker_filt3$absent_in_mutant[
  rownames(marker_filt3) %in% 
    c("npy2rl","fstl1b","slc1a2a","si:dkey-282h22.5")
] <- FALSE

marker_merge<-rbind(marker_filt,marker_filt3)
#write.csv(marker_merge,"c:/Users/admin/Desktop/Hb09_01_00 marker.csv")

#"filt---those are not tested in:
#"tac1","chrna4b", "desma","ENSDARG00000078898","lhfpl6","ENSDARG00000054537","igfbp5a", "mcm9","si:dkeyp-1h4.9","ENSDARG00000109844","si:dkey-183i3.6","atp2a2a","C17H14orf132","ENSDARG00000110426","ENSDARG00000096307","ENSDARG00000114719","ENSDARG00000088638","ahr1b","FNDC1","LY75"
#fiter3---those are not tested in :
#"zgc77752","fhl3b","ENSDARG00000067701","ENSDARG00000093936,","ENSDARG00000104624","fstl1b","si:dkey-238c7.12","BCO2","irf7","npy2rl","ENSDARG00000036462","ENSDARG00000100899","slc1a2a","si:dkey-282h22.5"

#saving marker genes in batch for mannual check
library(gridExtra)
for (x in 1:57){
  gene<-as.character(marker_merge[x,"gene"])
  p1<-FeaturePlot(WT,reduction = "umap_cca", gene)
  p2<-FeaturePlot(mut,reduction = "umap_cca", gene)
  p<-grid.arrange(p1,p2,ncol=2)
  dir<-file.path("C:/Users/admin/Desktop/Hb09_0_1 marker",paste0(gene,".png"))
  ggsave(dir,plot = p,width=16,height=8,dpi=500)
}

p1<-FeaturePlot(WT,reduction = "umap_cca", c("hnrnpm"))
p2<-FeaturePlot(mut,reduction = "umap_cca", c("hnrnpm"))

marker$gene<-rownames(marker)
for (x in 1:150){
  gene<-as.character(marker[x,"gene"])
  p1<-FeaturePlot(WT,reduction = "umap_cca", gene)
  p2<-FeaturePlot(mut,reduction = "umap_cca", gene)
  p<-grid.arrange(p1,p2,ncol=2)
  dir<-file.path("C:/Users/admin/Desktop/dont use",paste0(gene,".png"))
  ggsave(dir,plot = p,width=16,height=8,dpi=500)
}



