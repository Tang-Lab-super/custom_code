rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(harmony)

obj_WT_d14 <- readRDS("final.WT_d14.filtered.v2.rds")
obj_WT_d60 <- readRDS("final.WT_d60.filtered.v2.rds")
obj_KO_d14 <- readRDS("final.KO_d14.filtered.v2.rds")
obj_KO_d60 <- readRDS("final.KO_d60.filtered.v2.rds")

obj <- merge(obj_WT_d14, y = c(obj_WT_d30, obj_WT_d60, obj_KO_d14, obj_KO_d31, obj_KO_d60), merge.data = TRUE)

obj@meta.data$sample_type<-obj@meta.data$sample
obj@meta.data$sample_type<-gsub("_merge_5scRNA_.*","",obj@meta.data$sample_type)
obj@meta.data$sample_type<-gsub("Organoid_TLR_2","KO",obj@meta.data$sample_type)
obj@meta.data$sample_type<-gsub("Organoid_TLR","KO",obj@meta.data$sample_type)
obj@meta.data$sample_type<-gsub("Organoid_WT","WT",obj@meta.data$sample_type)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = row.names(obj))
obj <- RunPCA(object = obj, assay = "RNA", npcs = 20, features= VariableFeatures(obj)[1:2000])
theta1 = 0.5
pc = 20
obj <- RunHarmony(obj,group.by.vars=c("sample_type"),assay.use="RNA",reduction="pca",dims.use=1:pc,theta=c(theta1),
                  max.iter.harmony = 50, max.iter.cluster=200,kmeans_init_nstart=2, kmeans_init_iter_max=1000,
                  return_object = TRUE,plot_convergence = FALSE)
obj <- RunUMAP(obj, dims = 1:pc, reduction = "harmony", reduction.name = "umap_harmony", reduction.key = "UMAPharmony_")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:pc)
obj <- FindClusters(obj, verbose = T, resolution=0.6)

saveRDS(obj,"merge.harmony.rds")