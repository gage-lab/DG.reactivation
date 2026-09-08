library(Seurat)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(SeuratWrappers)
library(monocle3)

# Merge SMART-seq datasets
seurat.2026 <- read.csv("Data/SMARTseq_rawCounts_2026.csv",row.names=1)
pdata.2026 <- read.csv("Data/260310_SarahP_metadata_final.csv")
pdata.2026$sampleID <- paste0("X",pdata.2026$sampleID)
pdata.2026$batch <- "2026"

seurat.2025 <- read.csv("Data/SMARTseq_rawCounts_2025.csv",row.names=1)
pdata.2025 <- read.csv("Data/251121_SarahP_metadata_final.csv")
pdata.2025$sampleID <- paste0("X",pdata.2025$sampleID)
pdata.2025$batch <- factor(pdata.2025$sort_date)

seurat.2023 <- read.csv("Data/SMARTseq_rawCounts_2023.csv",row.names=1)
pdata.2023 <- read.csv("Data/230126_SarahP_metadata_final.csv")
pdata.2023$sampleID <- paste0("X",pdata.2023$sampleID)
pdata.2023$FOSraw <- NA
pdata.2023$GFPraw <- NA
pdata.2023$batch <- factor(pdata.2023$sort_plate)
levels(pdata.2023$batch) <- c("batch.1.2023","batch.1.2023","batch.1.2023","batch.1.2023","batch.2.2023","batch.2.2023","batch.2.2023","batch.1.2023","batch.2.2023")

# sort_plates (batch-effect):
# "35","36","37","76","Z7A" = batch.1
# "77","79","Z78","ZA1" = batch.2

df <- merge(seurat.2023,seurat.2025,by=0)
rownames(df) <- as.character(df$Row.names)
df <- df[,-1]

df <- merge(df,seurat.2026,by=0)
rownames(df) <- as.character(df$Row.names)
df <- df[,-1]

cols <- c("sampleID","mouse_number","mouse_birthdate","sex","exposure_type","timepoint","FullCategory","batch","FOSraw","GFPraw")

p <- rbind.data.frame(pdata.2026[,cols],pdata.2025[,cols],pdata.2023[,cols])
rownames(p) <- p$sampleID

p <- p[colnames(df),]

identical(rownames(p),colnames(df)) # TRUE

combined <- CreateSeuratObject(counts = df,min.cells = 1,meta.data=p) # 1920

# Preprocessing and normalization
combined <- combined %>%
  NormalizeData(verbose = FALSE)

all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

combined <- SplitObject(combined, split.by = "batch")

for (i in 1:length(combined)) {
  combined[[i]] <- SCTransform(combined[[i]], verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list = combined, nfeatures = 5000)
combined <- PrepSCTIntegration(object.list = combined, anchor.features = features)

combined <- FindIntegrationAnchors(object.list = combined, normalization.method = "SCT",anchor.features = features)
combined <- IntegrateData(anchorset = combined, normalization.method = "SCT")
DefaultAssay(combined) <- "integrated"

combined <- RunPCA(combined, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.15)

DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)


##### Figure S4A
DimPlot(combined)


##### Figure S4B
fivecolor <- c("blue","magenta2","green", "green4","orange2")
DimPlot(combined,group.by="FullCategory",cols=fivecolor)


##### Figure S4C
load("Data/DG_CA1_CA3_signature.rda") 

combined <- AddModuleScore(combined,features = list(ca1),name="CA1_sig")
combined <- AddModuleScore(combined,features = list(dg),name="DG_sig")
combined <- AddModuleScore(combined,features = list(ca3),name="CA3_sig")
FeaturePlot(combined,features = c("DG_sig1","CA1_sig1","CA3_sig1"),min.cutoff = "q5",order=T)

##### Figure S4D
VlnPlot(combined,features=c("nCount_RNA","nFeature_RNA"))


## Remove low quality cells and non-DG cells
# combined <- subset(combined, subset = seurat_clusters %in% c(0,1,2)) # 1688 cells
ids <- read.csv("Data/ids.csv",row.names=1)
combined <- subset(combined, cells = ids$x) # 1688 cells

## Redo clustering analysis
combined <- combined %>%
  NormalizeData(verbose = FALSE)

all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

combined <- SplitObject(combined, split.by = "batch")

for (i in 1:length(combined)) {
  combined[[i]] <- SCTransform(combined[[i]], verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list = combined, nfeatures = 5000)
combined <- PrepSCTIntegration(object.list = combined, anchor.features = features)

combined <- FindIntegrationAnchors(object.list = combined, normalization.method = "SCT",anchor.features = features)
combined <- IntegrateData(anchorset = combined, normalization.method = "SCT")
DefaultAssay(combined) <- "integrated"

combined <- RunPCA(combined, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.15)

DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)


##### Figure 3C
ncommscolor <- c("turquoise3", "chocolate1","gold1")
DimPlot(combined,cols=ncommscolor)


##### Figure 3D
DimPlot(combined,group.by="FullCategory",cols=fivecolor)


##### Figure 3E
p <- combined@meta.data
ggplot(p,aes(seurat_clusters,fill=FullCategory)) + 
    geom_bar(position = "fill") + 
    scale_fill_manual(values=fivecolor)

##### Figure 3F
ggplot(p,aes(seurat_clusters,fill=FullCategory)) + 
    geom_bar() + 
    scale_fill_manual(values=fivecolor)


##### Figure 3G
singleNE <- p[p$exposure_type %in% "singleNE",]
ggplot(singleNE,aes(FullCategory,fill=seurat_clusters)) + 
    geom_bar() + 
    scale_fill_manual(values=ncommscolor)


##### Figure 3H
reexp <- p[!p$exposure_type %in% "singleNE",]
ggplot(reexp,aes(FullCategory,fill=seurat_clusters)) + 
    geom_bar() + 
    scale_fill_manual(values=ncommscolor)

##### Figure 3I
combined@meta.data$exposure <- factor(combined@meta.data$exposure_type)
levels(combined@meta.data$exposure) <- c("reexposed","reexposed","reexposed","singleNE")

DimPlot(combined,group.by="exposure",cols=c('red','blue'))

p <- combined@meta.data
ggplot(p,aes(seurat_clusters,fill=exposure)) + 
    geom_bar(position="fill") + 
    scale_fill_manual(values=c('red','blue'))


##### Figure 3J
DimPlot(combined,group.by="timepoint",cols=c('purple','yellow','grey40'))    

p <- combined@meta.data
ggplot(p,aes(seurat_clusters,fill=timepoint)) + 
    geom_bar(position="fill") + 
    scale_fill_manual(values=c('purple','yellow','grey40'))

##### Figure 3K    
FeaturePlot(combined,features=c("Fos","Npas4","Arc","Homer1"),order=T)


##### Figure S4E
VlnPlot(combined,features=c("nCount_RNA","nFeature_RNA"))

##### Figure S4F    
FeaturePlot(combined,features=c("FOSraw"),order=T) + scale_color_gradientn(na.value = "white",colors=c("lightgrey", "blue"))
FeaturePlot(combined,features=c("GFPraw"),order=T) + scale_color_gradientn(na.value = "white",colors=c("lightgrey", "blue"))

##### Figure S4G  
tmp <- subset(combined, subset = !is.na(GFPraw))
VlnPlot(tmp,features=c("FOSraw","GFPraw"))


##### Differential expression analysis
reexp <- subset(combined, subset = exposure == "reexposed")
singleNE <- subset(combined, subset = exposure == "singleNE")

# Get the raw counts
cts <- GetAssayData(combined,assay = "RNA",layer = "counts")

### Comparison: Reactivated cluster 1 and compared GFP+ nuclei with any FOS status to GFP-FOS+ nuclei within Re-exposed mice

p <- reexp@meta.data 
p <- p[p$seurat_clusters %in% "1",]

p$groups <- factor(p$FullCategory)
levels(p$groups) <- c("GFPNegFOSPos","GFPPos","GFPPos","GFPPos")

df <- cts[,rownames(p)]
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = p,
                              design= ~ batch + groups)

dds <- DESeq(dds, test = "LRT", reduced = ~ batch, useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)  

res <- DESeq2::results(dds, contrast=c("groups","GFPPos","GFPNegFOSPos"))
res <- as.data.frame(res)

# No DEGs
sum(res$padj < 0.01 & !is.na(res$padj) & abs(res$log2FoldChange) > 1) 

##### Figure 4A
### Comparison: Recent Activation vs. Baseline [GFP-FOS+ in cluster 2 vs. GFP-FOS- in cluster 0 within Re-exposed]
p <- reexp@meta.data 
p <- p[(p$seurat_clusters %in% "2" & p$FullCategory %in% "GFPNegFOSPos") | (p$seurat_clusters %in% "0" & p$FullCategory %in% "GFPNegFOSNeg"),]

p$groups <- factor(p$FullCategory)

df <- cts[,rownames(p)]
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = p,
                              design= ~ batch + groups)

dds <- DESeq(dds, test = "LRT", reduced = ~ batch, useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)  

res <- DESeq2::results(dds, contrast=c("groups","GFPNegFOSPos","GFPNegFOSNeg"))
res <- as.data.frame(res)

ggplot(res,aes(x=log2FoldChange,y=-1*log10(padj),label=rownames(res))) + 
  geom_point() + 
  theme_classic() + 
  geom_text()

tmp <- combined@meta.data
tmp$groups <- "others"
tmp[rownames(p),"groups"] <- as.character(p$groups)
combined@meta.data <- tmp

DimPlot(combined,group.by="groups",cols=c('black','deeppink','grey'))


##### Figure 4C
### Comparison: Reactivation vs. Baseline [GFP+ or FOS+ in cluster 1 vs. GFP-FOS- in cluster 0 within Re-exposed]
p <- reexp@meta.data 
p <- p[(p$seurat_clusters %in% "1" & p$FullCategory %in% c("GFPPos","GFPPosFOSNeg","GFPPosFOSPos","GFPNegFOSPos")) | (p$seurat_clusters %in% "0" & p$FullCategory %in% "GFPNegFOSNeg"),]

p$groups <- factor(p$seurat_clusters)
p$groups <- relevel(p$groups,ref="0")

df <- cts[,rownames(p)]
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = p,
                              design= ~ batch + groups)

dds <- DESeq(dds, test = "LRT", reduced = ~ batch, useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)  

res <- DESeq2::results(dds, contrast=c("groups","1","0"))
res <- as.data.frame(res)

ggplot(res,aes(x=log2FoldChange,y=-1*log10(padj),label=rownames(res))) + 
  geom_point() + 
  theme_classic() + 
  geom_text()

tmp <- combined@meta.data
tmp$groups <- "others"
tmp[rownames(p),"groups"] <- as.character(p$groups)
combined@meta.data <- tmp

DimPlot(combined,group.by="groups",cols=c('black','goldenrod','grey'))


##### Figure 4E
### Comparison: Reactivation vs. Recent Activation [GFP+ in cluster 1 vs. GFP-FOS+ in cluster 2 within Re-exposed]
p <- reexp@meta.data 
p <- p[(p$seurat_clusters %in% "1" & p$FullCategory %in% c("GFPPos","GFPPosFOSNeg","GFPPosFOSPos")) | (p$seurat_clusters %in% "2" & p$FullCategory %in% "GFPNegFOSPos"),]

p$groups <- factor(p$seurat_clusters)
p$groups <- relevel(p$groups,ref="2")

df <- cts[,rownames(p)]
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = p,
                              design= ~ batch + groups)

dds <- DESeq(dds, test = "LRT", reduced = ~ batch, useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)  

res <- DESeq2::results(dds, contrast=c("groups","1","2"))
res <- as.data.frame(res)

ggplot(res,aes(x=log2FoldChange,y=-1*log10(padj),label=rownames(res))) + 
  geom_point() + 
  theme_classic() + 
  geom_text()

tmp <- combined@meta.data
tmp$groups <- "others"
tmp[rownames(p),"groups"] <- as.character(p$groups)
combined@meta.data <- tmp

DimPlot(combined,group.by="groups",cols=c('green','black','grey'))

##### Figure 4G
### Comparison: Remote Activation vs. Baseline [GFP+ in cluster 0 vs. GFP-FOS- in cluster 0 within Single NE]
p <- singleNE@meta.data 
p <- p[(p$seurat_clusters %in% "0" & p$FullCategory %in% c("GFPPos","GFPPosFOSNeg","GFPPosFOSPos")) | (p$seurat_clusters %in% "0" & p$FullCategory %in% "GFPNegFOSNeg"),]

p$groups <- factor(p$FullCategory)
levels(p$groups) <- c("GFPNegFOSNeg","GFPPos","GFPPos","GFPPos")

p$groups <- relevel(p$groups,ref="GFPNegFOSNeg")

df <- cts[,rownames(p)]
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = p,
                              design= ~ batch + groups)

dds <- DESeq(dds, test = "LRT", reduced = ~ batch, useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)  

res <- DESeq2::results(dds, contrast=c("groups","1","2"))
res <- as.data.frame(res)

ggplot(res,aes(x=log2FoldChange,y=-1*log10(padj),label=rownames(res))) + 
  geom_point() + 
  theme_classic() + 
  geom_text()

tmp <- combined@meta.data
tmp$groups <- "others"
tmp[rownames(p),"groups"] <- as.character(p$groups)
combined@meta.data <- tmp

DimPlot(combined,group.by="groups",cols=c('black','darkgreen','grey'))

##### Figure 4J
FeaturePlot(combined,features = c("Egr1","Nr4a2","Junb","Ptgs2"),min.cutoff = "q5",order=T)


##### Figure 4K
FeaturePlot(combined,features = c("Jarid2","Pam","Rnd3","Pdp1"),min.cutoff = "q5",order=T)


##### Figure 4L
FeaturePlot(combined,features = c("Penk","Acvr1c","Bhlhe41","Shisa4"),min.cutoff = "q5",order=T)

seurat.cds <- as.cell_data_set(combined)

seurat.cds <- cluster_cells(cds = seurat.cds)

seurat.cds <- learn_graph(seurat.cds)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(seurat.cds, type="GFPNegFOSNeg"){
  cell_ids <- which(colData(seurat.cds)[, "FullCategory"] == type)
  
  closest_vertex <-
    seurat.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(seurat.cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(seurat.cds)[["UMAP"]])$name[as.numeric(names
                                                                     (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
seurat.cds <- order_cells(seurat.cds, root_pr_nodes=get_earliest_principal_node(seurat.cds))

# Test all genes for spatial autocorrelation along the trajectory graph
pr_graph_test_res <- graph_test(seurat.cds, 
                                neighbor_graph = "principal_graph", 
                                cores = 10)

pr_graph_test_res <- pr_graph_test_res[with(pr_graph_test_res,order(q_value,-morans_I)),]

# Filter significant genes
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.01))

pseudotime <- data.frame(pseudotime=seurat.cds@principal_graph_aux@listData$UMAP$pseudotime)

identical(rownames(pseudotime),rownames(combined@meta.data))
combined$pseudotime <- pseudotime$pseudotime

##### Figure S5C
plot_cells(
  cds = seurat.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, cell_size = 1
)

combined@meta.data$groups <- factor(combined@meta.data$FullCategory)
combined@meta.data$groups <- factor(combined@meta.data$groups,levels(combined@meta.data$groups)[c(2,5,4,3,1)])

RidgePlot(combined,features = "pseudotime",group.by = "groups",cols = c("#e6194B","#f58231","#ffe119","#42d4f4","#4363d8"))

