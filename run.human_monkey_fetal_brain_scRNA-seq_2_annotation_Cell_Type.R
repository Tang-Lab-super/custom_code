# /home/yuchen/miniconda3/envs/R4.0/bin/R
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(harmony)
library(monocle3)
library(pheatmap)
library(RColorBrewer)
library(SeuratWrappers)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

obj <- readRDS("merge.harmony.rds")
pdf(file="DimPlot_by_cluster.pdf", height=10,width=11)
DimPlot(obj, pt.size=0.3, label=TRUE, label.size = 6)
dev.off()

markers <- FindMarkers(object = obj, ident.1=c(12), ident.2=c(2), slot='data', group.by="seurat_clusters", logfc.threshold=0.2, min.pct=0.1)
variable<- setdiff(VariableFeatures(obj), intersect(VariableFeatures(obj),row.names(head(markers[order(markers$p_val_adj,abs(markers$avg_log2FC)),],500))))
# variable<- VariableFeatures(obj)
# obj <- NormalizeData(obj)
# obj <- ScaleData(obj, features = row.names(obj))
obj <- RunPCA(object = obj, assay = "RNA", npcs = 30, features= variable)
obj <- RunHarmony(obj,group.by.vars=c("sample_type"),assay.use="RNA",reduction="pca",dims.use=1:30,theta=0.5,max.iter.harmony = 50, 
                  max.iter.cluster=200, kmeans_init_nstart=2, kmeans_init_iter_max=1000,return_object = TRUE,plot_convergence = FALSE)
obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony", reduction.key = "UMAPharmony_")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj)
pdf(file="0.DimPlot_by_cluster.1.pdf", height=10,width=11)
DimPlot(obj, pt.size=0.3, label=TRUE, label.size = 6)
dev.off()

obj <- subset(obj,sample_type %in% c("WT_d14","WT_d60","KO_d14","KO_d60"))
pdf(file="0.DimPlot_by_cluster.2.pdf", height=10,width=11)
DimPlot(obj, pt.size=0.3, label=TRUE, label.size = 6)
dev.off()

obj <- subset(obj,seurat_clusters %in% c(15,16,17,18,19,20,21,22),invert=T)
pdf(file="0.DimPlot_by_cluster.3.pdf", height=10,width=11)
DimPlot(obj, pt.size=0.3, label=TRUE, label.size = 6)
dev.off()

ko14<-readRDS("/data3/xiongdan/scRNA-seq/reRmDoublet/specificCells/final.KO_d14.filtered.v2.rds")
del<-sample(row.names(ko14@meta.data[ko14@meta.data$seurat_clusters %in% c(1),]),600)
del<-c(del,row.names(ko14@meta.data[ko14@meta.data$seurat_clusters %in% c(4,10),]))
saveRDS(del,"0.ko_d14.delete_cells.rds")
obj@meta.data$barcode<-row.names(obj@meta.data)
obj <- subset(obj,barcode %in% del,invert=T)
obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony", reduction.key = "UMAPharmony_")
obj <- FindNeighrm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(harmony)
library(monocle3)
library(pheatmap)
library(RColorBrewer)
library(SeuratWrappers)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

obj <- readRDS("merge.harmony.rds")
pdf(file="DimPlot_by_cluster.pdf", height=10,width=11)
DimPlot(obj, pt.size=0.3, label=TRUE, label.size = 6)
dev.off()

markers <- FindMarkers(object = obj, ident.1=c(12), ident.2=c(2), slot='data', group.by="seurat_clusters", logfc.threshold=0.2, min.pct=0.1)
variable<- setdiff(VariableFeatures(obj), intersect(VariableFeatures(obj),row.names(head(markers[order(markers$p_val_adj,abs(markers$avg_log2FC)),],500))))

obj <- RunPCA(object = obj, assay = "RNA", npcs = 30, features= variable)
obj <- RunHarmony(obj,group.by.vars=c("sample_type"),assay.use="RNA",reduction="pca",dims.use=1:30,theta=0.5,max.iter.harmony = 50, 
                  max.iter.cluster=200, kmeans_init_nstart=2, kmeans_init_iter_max=1000,return_object = TRUE,plot_convergence = FALSE)
obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony", reduction.key = "UMAPharmony_")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj)

obj <- subset(obj,sample_type %in% c("WT_d14","WT_d60","KO_d14","KO_d60"))
pdf(file="DimPlot_by_cluster.pdf", height=10,width=11)
DimPlot(obj, pt.size=0.3, label=TRUE, label.size = 6)
dev.off()


##################################################################################################################
obj@meta.data$type<-gsub("_d.*","",obj@meta.data$sample_type)
obj@meta.data$Cell_type<-""
obj@meta.data$Cell_type<-"unknow"
obj@meta.data[obj@meta.data$seurat_clusters %in% c(0,8,1,3,13,15),]$Cell_type<-"G2M"
obj@meta.data[obj@meta.data$seurat_clusters %in% c(2,6,11),]$Cell_type<-"NEC"
obj@meta.data[obj@meta.data$seurat_clusters %in% c(9,14),]$Cell_type<-"prog"
obj@meta.data[obj@meta.data$seurat_clusters %in% c(4),]$Cell_type<-"IPC"
obj@meta.data[obj@meta.data$seurat_clusters %in% c(7,10,12,16),]$Cell_type<-"EN"
obj@meta.data[obj@meta.data$seurat_clusters %in% c(5),]$Cell_type<-"IN"

pdf(file="DimPlot_by_celltype.pdf", height=10,width=11)
DimPlot(obj, group.by="Cell_type", pt.size=0.3, label=TRUE, label.size = 6)
dev.off()

obj@active.ident <- as.factor(obj$Cell_type)
cur_markers <- FindAllMarkers(obj, assays='RNA', slot='data', only.pos =TRUE, logfc.threshold =0.25)

library(plyr)
library(ggplot2)
data_celltype<-as.data.frame(table(obj@meta.data[,c("seurat_clusters","sample_type")]))
data_celltype<-merge(data_celltype,ddply(data_celltype,.(seurat_clusters),function(x){sum(x$Freq)}),by="seurat_clusters")
data_celltype$Frequency<-data_celltype$Freq/data_celltype$V1*100

data_celltype<-as.data.frame(table(obj@meta.data[,c("Cell_type","sample_type")]))
data_celltype<-merge(data_celltype,ddply(data_celltype,.(Cell_type),function(x){sum(x$Freq)}),by="Cell_type")
data_celltype$Frequency<-data_celltype$Freq/data_celltype$V1*100

data_celltype<-as.data.frame(table(obj@meta.data[,c("Cell_type","sample_type")]))
data_celltype<-merge(data_celltype,ddply(data_celltype,.(sample_type),function(x){sum(x$Freq)}),by="sample_type")
data_celltype$Frequency<-data_celltype$Freq/data_celltype$V1*100

saveRDS(obj,"harmony.obj.rds")

obj@reductions$umap<-obj@reductions$umap_harmony
cds <- as.cell_data_set(obj)
# clustering
cds <- cluster_cells(cds = cds, reduction_method = "UMAP", cluster_method="leiden", resolution=0.001)
cds <- learn_graph(cds = cds, use_partition = T, learn_graph_control=list(ncenter=1000))
# set start cell
startRoot <- sample(row.names(obj@meta.data[obj$Cell_type=="G2M",]),500)
cds <- order_cells(cds, reduction_method="UMAP", root_cells=startRoot)

saveRDS(cds,"monocle.cds.rds")
