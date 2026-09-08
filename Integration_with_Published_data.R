library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(DESeq2)

##### Integration with 4hr reactivation from Jaeger et al., 2018
# Reference: https://doi.org/10.1038/s41467-018-05418-8
# Raw counts and metadata were downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98679

load("Data/Jaeger_snRNAseq_rawCounts.rda")

jaeger.df <- df
jaeger.p <- pdata
jaeger.p <- jaeger.p[jaeger.p$cell.type %in% c("DG","Unknown"),]
jaeger.df <- jaeger.df[,rownames(jaeger.p)]
jaeger.p$batch <- "Jaeger"
jaeger.p$timepoint <- factor(jaeger.p$exposure)
levels(jaeger.p$timepoint) <- c("1hr","4hr","5hr","4hr","4hr","unk","unk","HC")
jaeger.p$exposure_type <- factor(jaeger.p$exposure)
jaeger.p$mouse_number <- factor(jaeger.p$strain)
jaeger.p$FullCategory <- factor(paste(jaeger.p$Fos,jaeger.p$Arc,sep="."))
levels(jaeger.p$FullCategory) <- c("FOSNegARCNeg","FOSNegARCPos","FOSPosARCNeg","FOSPosARCPos")
jaeger.p$cellID <- jaeger.p$sampleID

### Output from SMARTseq.R script
smartseq.counts <- as.data.frame(GetAssayData(combined,assay = "RNA",layer = "counts")) # raw counts
all <- merge(smartseq.counts,jaeger.df,by=0)
rownames(all) <- as.character(all$Row.names)
all <- all[,-1]

p <- combined@meta.data
p$study <- "Parylak"
jaeger.p$study <- "Jaeger"

cols <- c("cellID","batch","exposure_type","mouse_number","FullCategory","study")
pd <- rbind.data.frame(p[,cols],jaeger.p[,cols])

seurat <- CreateSeuratObject(counts = all,min.cells = 1,meta.data=pd) # 2323 cells

seurat <- seurat %>%
  NormalizeData(verbose = FALSE)

all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

seurat <- SplitObject(seurat, split.by = "batch")

for (i in 1:length(seurat)) {
  seurat[[i]] <- SCTransform(seurat[[i]], verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list = seurat, nfeatures = 5000)
seurat <- PrepSCTIntegration(object.list = seurat, anchor.features = features)

seurat <- FindIntegrationAnchors(object.list = seurat, normalization.method = "SCT",anchor.features = features)
seurat <- IntegrateData(anchorset = seurat, normalization.method = "SCT")
DefaultAssay(seurat) <- "integrated"

seurat <- RunPCA(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30, verbose = FALSE)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
seurat <- FindClusters(seurat, resolution = 0.15)

DefaultAssay(seurat) <- "RNA"
seurat <- JoinLayers(seurat)

##### Figure S5D
DimPlot(seurat)

VlnPlot(seurat,features = "Penk",split.by = "study")

##### Figure S5E
FeaturePlot(seurat,features = c("Sorcs3","Penk","Arc","Fos"),order=T,min.cutoff = "q5")

##### Integration with Rao-Ruiz et al. PMID:31110186
## Data: Table S3, Figure 2C

g <- read.csv("Data/Rao.Ruiz.Table_S3.csv")
genes <- g[g$log2FoldChange > 0,"Gene.name"]

# alternative gene symbols
alternative <- read.table(text="Lphn3   Adgrl3  adhesion G protein-coupled receptor L3  MGI:2441950
AI414108    Igsf9b  immunoglobulin superfamily, member 9B   MGI:2685354
1190002N15Rik   Dipk2a  divergent protein kinase domain 2A  MGI:1916111
Rfwd2   Cop1    COP1, E3 ubiquitin ligase   MGI:1347046
Ppap2a  Plpp1   phospholipid phosphatase 1  MGI:108412
Epb4.1l2    Epb41l2 erythrocyte membrane protein band 4.1 like 2    MGI:103009
Tmem2   Cemip2  cell migration inducing hyaluronidase 2 MGI:1890373
Lepre1  P3h1    prolyl 3-hydroxylase 1  MGI:1888921
2610301B20Rik   Cfap418 cilia and flagella associated protein 418   MGI:1914407
Lins    Lins1   lines homolog 1 MGI:1919885
Grasp   Tamalin trafficking regulator and scaffold protein tamalin  MGI:1860303
Pqlc2   Slc66a1 solute carrier family 66 member 1   MGI:2384837
Ctage5  Mia2    MIA SH3 domain ER export factor 2   MGI:2159614
Prosc   Plpbp   pyridoxal phosphate binding protein MGI:1891207
Mki67ip Nifk    nucleolar protein interacting with the FHA domain of MKI67  MGI:1915199
Vimp    Selenos selenoprotein S MGI:95994
2610019F03Rik   Tdrp    testis development related protein  MGI:1919398
Tmem48  Ndc1    NDC1 transmembrane nucleoporin  MGI:1920037
Ccdc90a Mcur1   mitochondrial calcium uniporter regulator 1 MGI:1923387
Ppapdc1a    Plpp4   phospholipid phosphatase 4  MGI:2685936
AI314180    Ecpas   Ecm29 proteasome adaptor and scaffold   MGI:2140220
Ufd1l   Ufd1    ubiquitin recognition factor in ER-associated degradation 1 MGI:109353",sep="\t")
genes <- c(genes,alternative$V2)

combined <- AddModuleScore(combined,features = list(genes),name = "Rao_Ruiz_Table_S3.")

##### Figure S5E
tmp <- combined@meta.data
tmp$groups <- as.character(tmp$seurat_clusters)
tmp[tmp$seurat_clusters %in% "0","groups"] <- "Baseline, GFP+"
tmp[tmp$seurat_clusters %in% "0" & tmp$FullCategory %in% "GFPNegFOSNeg","groups"] <- "Baseline, GFP-"
tmp$groups <- factor(tmp$groups)
levels(tmp$groups)[1:2] <- c("Reactivated","Newly active")
tmp$groups <- factor(tmp$groups, levels=c("Baseline, GFP-","Baseline, GFP+","Newly active","Reactivated"))
combined@meta.data <- tmp

FeaturePlot(combined,features = "Rao_Ruiz_Table_S3.1",order=T,min.cutoff = "q10")

VlnPlot(combined,features = "Rao_Ruiz_Table_S3.1",group.by = "groups")


##### Figure S5E
FeaturePlot(seurat,features = c("Sorcs3","Penk","Arc","Fos"),order=T,min.cutoff = "q5")

##### Integration with Chen et al. PMID: 33177708
## Data: Figure 3C

g <- read.table(text="Hid1
Kctd10
Pdha1
Pigq
Mtrex
Abcf3
Miga2
Mmd
Sar1a
Eif2ak1
Acsf3
Cdc42se2
Retreg2
Vamp2
Gsk3b
Hnrnph2
Pja2
Sdha
Gdi2
Rab15
Fam131a
Gfra2
Pip4k2c
Ncdn
Usp5
Nell2
Tmem151a
Dmtn
Rtn3
Pcsk2
Pfkm
Trim32
Mfsd14b
Rab5a
Tdg
Emc1
Gpm6a
Elmod1
Cycs
Sarnp
Rab24
Erp29
Ghitm
Zfp706
Sult4a1
App
Cck
Emc4
Psmb6
Atp6v0b
Hint1
Guk1
Rtn1
Clpp
Crip2
Mpc1
Atp6v0c
Atp5g3
Tmem50a
Plekhb2
Syt13
Garnl3
Dpysl4
Aplp1
Hnrnpk
Nsf
Mfsd14a
Nck2
Stx1b
Pak1
Slc25a46
Itfg1
Lmbrd1
Tmx1
Dner
Atad1
Ankrd45
Timm29
Vopp1
Pls3
Hmg20a
Ctbp1
Strip1
Cdv3
Inpp5f
Prkar1b
Slc30a9
Alg2
Trim35
Hacd3
Serinc1
Serinc3
Ptp4a1")

combined <- AddModuleScore(combined,features = list(g$V1),name = "Chen_Fig3c.")

##### Figure S5F

FeaturePlot(combined,features = "Chen_Fig3c.1",order=T,min.cutoff = "q10")

VlnPlot(combined,features = "Chen_Fig3c.1",group.by = "groups")

##### Integration with Marco et al. PMID 33020654
## Data: Figure 5A

g <- read.csv("Data/Marco_SuppTable8_Late_vs_Reactivated.csv")
genes <- g[g$log2FoldChange > 0,"ensembl_gene_id"]
alternative <- read.table(text="Hist1h2bc   H2bc4   H2B clustered histone 4 MGI:1915274
Hist1h2af   H2ac10  H2A clustered histone 10    MGI:2448309
Hist1h4h    H4c8    H4 clustered histone 8  MGI:2448427
Gpr98   Adgrv1  adhesion G protein-coupled receptor V1  MGI:1274784
Usmg5   Atp5mk  ATP synthase membrane subunit k MGI:1891435
Cxx1a   Rtl8a   retrotransposon Gag like 8A MGI:1913408
Ccdc135 Drc7    dynein regulatory complex subunit 7 MGI:2685616
Gm1821  Ubb-ps  ubiquitin B, pseudogene MGI:3037679
2900055J20Rik   Kctd16  potassium channel tetramerisation domain containing 16  MGI:1914659
2610019F03Rik   Tdrp    —   —
Ppapdc1a    Plpp4   phospholipid phosphatase 4  MGI:2685936
LOC102634431    Gm32014 predicted gene, 32014   —",sep="\t")
genes <- c(genes,alternative$V2)

combined <- AddModuleScore(combined,features = list(genes),name = "Marco_Table_S8.Late.vs.Reactivated.")

##### Figure S5G

FeaturePlot(combined,features = "Marco_Table_S8.Late.vs.Reactivated.1",order=T,min.cutoff = "q10")

VlnPlot(combined,features = "Marco_Table_S8.Late.vs.Reactivated.1",group.by = "groups")
