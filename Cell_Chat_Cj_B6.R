library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
library(patchwork)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)
B6 <- readRDS("H:/Manuscripts/Mucosal Immunology Paper/Revised Figures/Cj_B6_Seurat_Analysis.rds")
DimPlot(B6, reduction = "umap", label = TRUE, pt.size = 0.5)
data.input <- GetAssayData(B6,assay = "RNA",slot = "counts")
B6_mod <- CreateSeuratObject(counts = data.input)
B6_mod <- NormalizeData(B6_mod)
levels(B6@meta.data$seurat_clusters)
Idents(B6) <- 'seurat_clusters'
new.seurat_clusters<- c("MP_MHCII-","Neut.LCN2-low", "Neut.Lcn2-hi", "Mono",  "MP", "NK","T_cells", "GD_T", "ILC2", "Mast cells", "cDC1", "Cluster-11", "B", "mDC") 
names(new.seurat_clusters) <- levels(B6)
B6
B6 <- RenameIdents(B6, new.seurat_clusters)
B6_mod$celltypes <- B6@active.ident
meta<-B6_mod@meta.data
cellchat <- createCellChat(object = data.input, meta=meta, group.by = "celltypes")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
saveRDS(cellchat, file='H:/Manuscripts/Mucosal Immunology Paper/Revised Figures/B6_cellchat.rds')
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
cellchat@netP$pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
