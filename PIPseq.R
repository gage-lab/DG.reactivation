library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(DESeq2)

pipseq.counts <- read.csv("Data/PIPseq_rawCounts.csv",row.names = 1)
pdata <- read.csv("Data/PIPseq_metadata.csv",row.names = 1)

fourcolor <- c("deeppink","green4","black","#e08c05")

pipseq <- CreateSeuratObject(counts = pipseq.counts,min.cells = 0,meta.data=pdata,min.features = 0) # 5887 cells

pipseq <- NormalizeData(pipseq, normalization.method = "LogNormalize", scale.factor = 10000)
pipseq <- FindVariableFeatures(pipseq, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pipseq)
pipseq <- ScaleData(pipseq, features = all.genes)

##Run PCA
pipseq <- RunPCA(pipseq, features = VariableFeatures(object = pipseq))

pipseq <- FindNeighbors(pipseq, dims = 1:25)
pipseq <- RunUMAP(pipseq, dims = 1:25)
pipseq <- FindClusters(pipseq, resolution = 0.25)

##### Figure S5H
DimPlot(pipseq,group.by = "orig.ident")

##### Figure S5I
VlnPlot(pipseq,features = c("nFeature_RNA","nCount_RNA"))

##### Figure 5C
DimPlot(pipseq,group.by = "seurat_clusters",label = T,label.box = T,cols=fourcolor)

##### Figure 5D
FeaturePlot(pipseq,features = c("Fos","Npas4","Arc","Homer1","Penk","Pdp1","Shisa4","Itgav","Acvr1c","Sorcs3","Rgs6","Rasgrp1"),ncol = 4,order=T,min.cutoff = "q5")

##### Figure 5E
DimPlot(pipseq,group.by = "exposure",cols=c("red","blue"))

##### Figure 5F
p <- pipseq@meta.data

PenkClass <- WhichCells(pipseq, expression = Penk > 0)
p$PenkClass <- "Negative"
p[PenkClass,"PenkClass"] <- "Positive"

ggplot(p,aes(PenkClass,fill=seurat_clusters)) + 
  geom_bar(position="fill") +
  scale_fill_manual(values=fourcolor,name="Cluster") + 
  labs(x="Penk Expression",y="Percentage (%)")

ggplot(p,aes(PenkClass,fill=seurat_clusters)) + 
  geom_bar() +
  scale_fill_manual(values=fourcolor,name="Cluster") + 
  labs(x="Penk Expression",y="Number of cells")

### Create pseudobulk profiles for differential analysis
pseudo_PIPseq <- AggregateExpression(pipseq, assays = "RNA", return.seurat = T, group.by = c("orig.ident", "seurat_clusters"))
Idents(pseudo_PIPseq) <- "seurat_clusters"

##### Figure 5I
### Recently active vs. Baseline
pseudodegs.recent.baseline.all <- FindMarkers(object = pseudo_PIPseq, 
                                              ident.1 = "0", 
                                              ident.2 = c("1","2"),
                                              test.use = "DESeq2")

res <- na.omit(pseudodegs.recent.baseline.all)
res$direction <- "NS"
res[res$avg_log2FC >= 1 & res$p_val_adj < 0.01, "direction"] <- "UP"
res[res$avg_log2FC <= -1 & res$p_val_adj < 0.01, "direction"] <- "DOWN"
res$gene <- rownames(res)
Top_Hits = head(arrange(res,p_val_adj),20)$gene 

ggplot(res, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) + 
  geom_hline(yintercept = -log10(0.01),  linetype = "dashed", col = "gray70") + geom_vline(xintercept = c(-1,1), linetype = "dashed", col = "gray70") + 
  geom_point() + scale_color_manual(values = c("gray20","gray70","deeppink"), name = "") + 
  coord_cartesian(xlim = c(-10,10)) + 
  geom_text_repel(data=res[res$gene %in% Top_Hits,],aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction,label=gene)) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12))
  
##### Figure 5J
### Reactivated vs. Baseline
pseudodegs.reactivated.baseline.all <- FindMarkers(object = pseudo_PIPseq, 
                                                   ident.1 = "3", 
                                                   ident.2 = c("1","2"),
                                                   test.use = "DESeq2")

res <- na.omit(pseudodegs.reactivated.baseline.all)
res$direction <- "NS"
res[res$avg_log2FC >= 1 & res$p_val_adj < 0.01, "direction"] <- "UP"
res[res$avg_log2FC <= -1 & res$p_val_adj < 0.01, "direction"] <- "DOWN"
res$gene <- rownames(res)
Top_Hits = head(arrange(res,p_val_adj),20)$gene 

ggplot(res, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) + 
  geom_hline(yintercept = -log10(0.01),  linetype = "dashed", col = "gray70") + geom_vline(xintercept = c(-1,1), linetype = "dashed", col = "gray70") + 
  geom_point() + scale_color_manual(values = c("gray20","gray70","#e08c05"), name = "") + 
  coord_cartesian(xlim = c(-10,10)) + 
  geom_text_repel(data=res[res$gene %in% Top_Hits,],aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction,label=gene)) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12))  
  
##### Figure 5K
### Reactivated vs. Recent
pseudodegs.reactivated.recent.all <- FindMarkers(object = pseudo_PIPseq, 
                                                 ident.1 = "3", 
                                                 ident.2 = "0",
                                                 test.use = "DESeq2")

res <- na.omit(pseudodegs.reactivated.recent.all)
res$direction <- "NS"
res[res$avg_log2FC >= 1 & res$p_val_adj < 0.01, "direction"] <- "UP"
res[res$avg_log2FC <= -1 & res$p_val_adj < 0.01, "direction"] <- "DOWN"
res$gene <- rownames(res)
Top_Hits = head(arrange(res,p_val_adj),20)$gene 

ggplot(res, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) + 
  geom_hline(yintercept = -log10(0.01),  linetype = "dashed", col = "gray70") + geom_vline(xintercept = c(-1,1), linetype = "dashed", col = "gray70") + 
  geom_point() + scale_color_manual(values = c("gray20","gray70","green3"), name = "") + 
  coord_cartesian(xlim = c(-10,10)) + 
  geom_text_repel(data=res[res$gene %in% Top_Hits,],aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction,label=gene)) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12))  

##### Figure 5L
### Baseline1 vs Baseline2
pseudodegs.baseline1.baseline2.all <- FindMarkers(object = pseudo_PIPseq, 
                                                  ident.1 = "1", 
                                                  ident.2 = "2",
                                                  test.use = "DESeq2")
res <- na.omit(pseudodegs.baseline1.baseline2.all)
res$direction <- "NS"
res[res$avg_log2FC >= 1 & res$p_val_adj < 0.01, "direction"] <- "UP"
res[res$avg_log2FC <= -1 & res$p_val_adj < 0.01, "direction"] <- "DOWN"
res$gene <- rownames(res)
Top_Hits = head(arrange(res,p_val_adj),20)$gene 

ggplot(res, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) + 
  geom_hline(yintercept = -log10(0.01),  linetype = "dashed", col = "gray70") + geom_vline(xintercept = c(-1,1), linetype = "dashed", col = "gray70") + 
  geom_point() + scale_color_manual(values = c("gray20","gray70","green4"), name = "") + 
  coord_cartesian(xlim = c(-10,10)) + 
  geom_text_repel(data=res[res$gene %in% Top_Hits,],aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction,label=gene)) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12))  
