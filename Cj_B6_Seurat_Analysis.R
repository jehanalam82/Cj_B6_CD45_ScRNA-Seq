remotes::install_version("Seurat", version = "3.X.X")
source("https://z.umn.edu/archived-seurat")
remotes::install_version("Seurat", version = "3.X.X")
install.packages('Seurat', version = "3.X.X")
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
B6

B6.data <- Read10X(data.dir = "H:/Manuscripts/ssRNA_Seq_B6_Normal_Conjuctiva/Single cell Analysis/NS-B6-18-20WEEKS/B6")
B6 <- CreateSeuratObject(counts = B6.data, project = "B6", min.cells = 3, min.features = 200)
B6
B6[["percent.mt"]] <- PercentageFeatureSet(B6, pattern = "^mt-")
#####Fig. S1A##########
VlnPlot(B6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(B6, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(B6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
B6 <- subset(B6, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
B6
B6 <- NormalizeData(B6, normalization.method = "LogNormalize", scale.factor = 10000)
B6 <- NormalizeData(B6)
B6 <- FindVariableFeatures(B6, selection.method = "vst", nfeatures = 2500)
top15 <- head(VariableFeatures(B6), 15)
plot1 <- VariableFeaturePlot(B6)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
plot1 + plot2
all.genes <- rownames(B6)
B6 <- ScaleData(B6, features = all.genes)
B6 <- RunPCA(B6, features = VariableFeatures(object = B6))
print(B6[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(B6, dims = 1:2, reduction = "pca")
B6 <- JackStraw(B6, num.replicate = 100)
B6 <- ScoreJackStraw(B6, dims = 1:20)
JackStrawPlot(B6, dims = 1:15)
ElbowPlot(B6)
B6 <- FindNeighbors(B6, dims = 1:20)
B6 <- FindClusters(B6, resolution = 0.5)
head(Idents(B6), 5)
B6 <- RunUMAP(B6, dims = 1:20)
DimPlot(B6, reduction = "umap")
#####################Fig.S1B#################
FeaturePlot(B6, features = c("Serpinb2","Lcn2","S100a9","Apoe","Cd209a", "Gzma","Cd28","Trdc","Il5","Xcr1","Mcpt4", "Cd79a", "Fscn1","Retnla", "Cd8a", "Cd4", "Foxp3"))
################Fig. 1B################
new.cluster.ids <- c("MP_MHCII-","Neut.LCN2-low", "Neut.Lcn2-hi", "Mono",  "MP", "NK","T_cells", "GD_T", "ILC2", "Mast cells", "cDC1", "Cluster-11", "B", "mDC")
names(new.cluster.ids) <- levels(B6)
B6 <- RenameIdents(B6, new.cluster.ids)
DimPlot(B6, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() 
DefaultAssay(B6) <- "RNA"
########cells cluster count#########
B6[["my.clusters"]] <- Idents(object = B6)
table(object@meta.data$my.clusters, object@meta.data$orig.ident)
###############S.Table.2#################
B6.markers <- FindAllMarkers(B6, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(B6.markers, file = "H:/Manuscripts/Mucosal Immunology Paper/Revised Figures/B6_Revision.markers" )
#################Fig.S2############
B6.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(B6, features = top10$gene, angle = 90 )
#####Fig.2A
features <- c("Itgam","Itgax","Ccr2")
RidgePlot(B6, features = features, ncol = 3)
######Fig. 2B
features <- c("Serpinb2","F5","Alox15","Ccrl2","Nlrp3","Marcksl1","S100a9", "S100a8","Lrg1","Apoe","Lyz2","Mafb","Cd209a","Cd74","H2-Ab1", "Gzma","Ccl5","AW112010","Ctla4","Cd28","Icos","Il17a","Trdc","Cxcr6","Il5","Gata3","Areg", "Mcpt4","Cpa3","Cma1","Ifi205","Cst3","Naaa","Jund","Egr1","Fos","Igkc","Cd79a","Ebf1","Ccl22","Ccr7","Tbc1d4"  )
DotPlot(B6, features = features, cols = c("blue", "red"), dot.scale = 6) + RotatedAxis()
##########Fig.4, Fig.3D and Fig.S3#######
VlnPlot(B6, features = c("Tnf"))+ NoLegend()
####"Fscn1", "Ccr7", "Cst3","H2-Ab1","Apoe","Csf1r","Cd209a","Serpinb2","Adgre1","Il5","Il13","Rxra","Aldh1a2","Thbs1","Arg1","Il13","Il10","Il1a","Il1b"
####"Il12a", "Il12b","Il18","Il6","Nos2","Ifng","Cd83","H2-Ab1","Casp1","Nlrp3","S100a9","Mmp9","Il17a","Tnf"
#########################################################
saveRDS(B6, file = "H:/Manuscripts/Mucosal Immunology Paper/Revised Figures/Cj-B6-NS.rds")
