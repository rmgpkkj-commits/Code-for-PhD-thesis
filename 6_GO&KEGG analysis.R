#- Project 'E:/seurat/fgsea' loaded. [renv 1.1.4]
library(enrichplot)
library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(org.Dr.eg.db)
library(clusterProfiler)

#=============================GSEA======================================================
data <- read_csv("C:/Users/admin/Desktop/2025-11-25 correction/DEseq result/Hb15_WT_Deseq.csv")  # Replace with your actual file name
# Make sure gene symbols are character
data$gene <- as.character(data$...1)

# rank by logFC
# Named vector of logFC (names = gene symbols)
geneList <- data$log2FoldChange   #(negative )
names(geneList) <- data$gene
# Sort by decreasing logFC
geneList <- sort(geneList, decreasing = TRUE) 

# For GSEA
geneList_entrez <- bitr(names(geneList), fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Dr.eg.db)
geneList_named <- geneList[geneList_entrez$SYMBOL]
names(geneList_named) <- geneList_entrez$ENTREZID

gse <- gseGO(geneList = geneList,    
             OrgDb = org.Dr.eg.db,
             keyType = "SYMBOL",
             ont = "ALL",
             minGSSize = 10,
             maxGSSize = 500,
             pvalueCutoff = 0.5)

gseKEGG_res <- gseKEGG(
  geneList     = geneList_named,   # must be ENTREZID-named vector
  organism     = "dre",            # zebrafish KEGG code
  keyType      = "ncbi-geneid",    # KEGG uses ENTREZ
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.5,
  verbose      = FALSE
)
# Convert Entrez IDs in leadingEdge to gene symbols
gseKEGG_res@result$core_enrichment <- sapply(gseKEGG_res@result$core_enrichment, function(x) {
  entrez <- unlist(strsplit(x, "/"))
  symbols <- mapIds(org.Dr.eg.db, keys = entrez, column = "SYMBOL",
                    keytype = "ENTREZID", multiVals = "first")
  paste(symbols, collapse = "/")
})

go_Hb15_full<-gse
kegg_Hb15_full<-gseKEGG_res

#============================summarizing and saveing results==================================================
#for cluster specific comparsion
GO_cluster_fulldata<-list(go_Hb02_full,go_Hb04_full,go_Hb05_full,go_Hb06_full,go_Hb07_full,
                          go_Hb08_full,go_Hb10_full,go_Hb15_full)
names(GO_cluster_fulldata) <- c(
  "Hb02","Hb04","Hb05","Hb06","Hb07","Hb08","Hb10","Hb15"
)
KEGG_cluster_fulldata<-list(kegg_Hb02_full,kegg_Hb04_full,kegg_Hb05_full,kegg_Hb06_full,kegg_Hb07_full,
                            kegg_Hb08_full,kegg_Hb10_full,kegg_Hb15_full)
names(KEGG_cluster_fulldata) <- c(
  "Hb02","Hb04","Hb05","Hb06","Hb07","Hb08","Hb10","Hb15"
)
#saveRDS(GO_cluster_fulldata,file = "C:/Users/admin/Desktop/GO_cluster_fulldata.rds")
#saveRDS(KEGG_cluster_fulldata,file = "C:/Users/admin/Desktop/KEGG_cluster_fulldata.rds")


#for cluster marker
GO_cluster_fulldata<-list(go_Hb01_full,go_Hb02_full,go_Hb03_full,go_Hb04_full,go_Hb05_full,go_Hb06_full,go_Hb07_full,
                          go_Hb08_full,go_Hb09_full,go_Hb10_full,go_Hb11_full,go_Hb13_full,go_Hb15_full)
names(GO_cluster_fulldata) <- c(
  "Hb01","Hb02","Hb03","Hb04","Hb05","Hb06","Hb07","Hb08","Hb09","Hb10","Hb11","Hb13","Hb15"
)
KEGG_cluster_fulldata<-list(kegg_Hb01_full,kegg_Hb02_full,kegg_Hb03_full,kegg_Hb04_full,kegg_Hb05_full,kegg_Hb06_full,kegg_Hb07_full,
                            kegg_Hb08_full,kegg_Hb09_full,kegg_Hb10_full,kegg_Hb11_full,kegg_Hb13_full,kegg_Hb15_full)
names(KEGG_cluster_fulldata) <- c(
  "Hb01","Hb02","Hb03","Hb04","Hb05","Hb06","Hb07","Hb08","Hb09","Hb10","Hb11","Hb13","Hb15"
)
saveRDS(GO_cluster_fulldata,file = "C:/Users/admin/Desktop/GO_clustermarker_fulldata.rds")
saveRDS(KEGG_cluster_fulldata,file = "C:/Users/admin/Desktop/KEGG_clustermarker_fulldata.rds")

#===========================summarizing and saveing results(continue)====================================================
# for per cluster comparsion between genotype
go_Hb02 <- go_Hb02_full@result
go_Hb04 <- go_Hb04_full@result
go_Hb05 <- go_Hb05_full@result
go_Hb06 <- go_Hb06_full@result
go_Hb07 <- go_Hb07_full@result
go_Hb08 <- go_Hb08_full@result
go_Hb10 <- go_Hb10_full@result
go_Hb15 <- go_Hb15_full@result
kegg_Hb02 <- kegg_Hb02_full@result
kegg_Hb04 <- kegg_Hb04_full@result
kegg_Hb05 <- kegg_Hb05_full@result
kegg_Hb06 <- kegg_Hb06_full@result
kegg_Hb07 <- kegg_Hb07_full@result
kegg_Hb08 <- kegg_Hb08_full@result
kegg_Hb10 <- kegg_Hb10_full@result
kegg_Hb15 <- kegg_Hb15_full@result

list_go <- list(
  go_Hb02, go_Hb04, go_Hb05,
  go_Hb07, go_Hb08, go_Hb10,go_Hb15)
cluster_names <- c("Hb02","Hb04","Hb05","Hb07","Hb08","Hb10","Hb15")
# Add cluster names
for (i in seq_along(list_go)) {
  list_go[[i]]$cluster <- cluster_names[i]
}
# Combine into one dataframe
go <- do.call(rbind, list_go)

list_kegg <- list(kegg_Hb04,kegg_Hb05,kegg_Hb06,kegg_Hb07,kegg_Hb08,
                  kegg_Hb10,kegg_Hb15)
cluster_names <- c("kegg_Hb04","kegg_Hb05","kegg_Hb06","kegg_Hb07","kegg_Hb08","kegg_Hb10","kegg_hb15")
# Add cluster names
for (i in seq_along(list_kegg)) {
  list_kegg[[i]]$cluster <- cluster_names[i]
}
# Combine into one dataframe
kegg <- do.call(rbind, list_kegg)

#write.csv(go,file = "C:/Users/admin/Desktop/GO_cluster.csv")
#write.csv(kegg,file = "C:/Users/admin/Desktop/kegg_cluster.csv")

#for cluster marker 
go_Hb01 <- go_Hb01_full@result
go_Hb02 <- go_Hb02_full@result
go_Hb03 <- go_Hb03_full@result
go_Hb04 <- go_Hb04_full@result
go_Hb05 <- go_Hb05_full@result
go_Hb06 <- go_Hb06_full@result
go_Hb07 <- go_Hb07_full@result
go_Hb08 <- go_Hb08_full@result
go_Hb09 <- go_Hb09_full@result
go_Hb10 <- go_Hb10_full@result
go_Hb11 <- go_Hb11_full@result
go_Hb13 <- go_Hb13_full@result
go_Hb15 <- go_Hb15_full@result
kegg_Hb01 <- kegg_Hb01_full@result
kegg_Hb02 <- kegg_Hb02_full@result
kegg_Hb03 <- kegg_Hb03_full@result
kegg_Hb04 <- kegg_Hb04_full@result
kegg_Hb05 <- kegg_Hb05_full@result
kegg_Hb06 <- kegg_Hb06_full@result
kegg_Hb07 <- kegg_Hb07_full@result
kegg_Hb08 <- kegg_Hb08_full@result
kegg_Hb09 <- kegg_Hb09_full@result
kegg_Hb10 <- kegg_Hb10_full@result
kegg_Hb11 <- kegg_Hb11_full@result
kegg_Hb13 <- kegg_Hb13_full@result
kegg_Hb15 <- kegg_Hb15_full@result

list_go <- list(
  go_Hb01,go_Hb02,go_Hb03, go_Hb04, go_Hb05,go_Hb06,
  go_Hb07, go_Hb08, go_Hb09,go_Hb10, go_Hb11,go_Hb13, go_Hb15
)
cluster_names <- c("Hb01","Hb02","Hb03","Hb04","Hb05","Hb06","Hb07","Hb08","Hb09","Hb10","Hb11","Hb13","Hb15")
for (i in seq_along(list_go)) {
  list_go[[i]]$cluster <- cluster_names[i]
}
go <- do.call(rbind, list_go)

list_kegg <- list(kegg_Hb01,kegg_Hb02,kegg_Hb03,kegg_Hb04,kegg_Hb05,kegg_Hb06,kegg_Hb07,kegg_Hb08,
                  kegg_Hb09,kegg_Hb10,kegg_Hb11,kegg_Hb13,kegg_Hb15)
cluster_names <- c("Hb01","Hb02","Hb03","Hb04","Hb05","Hb06","Hb07","Hb08","Hb09","Hb10","Hb11","Hb13","Hb15")
for (i in seq_along(list_kegg)) {
  list_kegg[[i]]$cluster <- cluster_names[i]
}
kegg <- do.call(rbind, list_kegg)

#write.csv(go,file = "C:/Users/admin/Desktop/GO_clustermarker.csv")
#write.csv(kegg,file = "C:/Users/admin/Desktop/kegg_clustermarker.csv")

#===============================annotations================================================
go_clust<-read.csv("C:/Users/admin/Desktop/GO_clustermarker.csv")
kegg_clust<-read.csv("C:/Users/admin/Desktop/kegg_clustermarker.csv")
go<-read.csv("C:/Users/admin/Desktop/GO_cluster.csv")
kegg<-read.csv("C:/Users/admin/Desktop/kegg_cluster.csv")

go<-go[go$p.adjust<0.05,]
kegg<-kegg[kegg$p.adjust<0.05,]

go_clust<-go_clust[go_clust$p.adjust<0.05,]
kegg_clust<-kegg_clust[kegg_clust$p.adjust<0.05,]



