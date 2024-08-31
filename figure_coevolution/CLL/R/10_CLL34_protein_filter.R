# filter Seurat object

library(dplyr)
library(ggplot2)
library(Seurat)

so.12 <- readRDS("./data/coevolution/objects/20220504_CLL34.rds")

# remove doublets and unannotated cells
annotation <- read.csv2("./data/coevolution/CLL3/20220506_deconvolution.csv", header = T, row.names = 1) %>% as.data.frame()
annotation <- rbind(annotation, read.csv2("./data/coevolution/CLL4/20220506_deconvolution.csv", header = T, row.names = 1) %>% as.data.frame())

so.12 <- subset(so.12, cells = annotation$bc[which(!is.na(annotation$annotation))])

so.12 <- NormalizeData(so.12, normalization.method = "CLR", margin = 2)
so.12 <- FindVariableFeatures(so.12)
so.12 <- ScaleData(so.12, features = rownames(so.12))
so.12 <- RunPCA(so.12)
so.12 <- RunUMAP(so.12, dims = 1:15)
so.12 <- FindNeighbors(so.12)
so.12 <- FindClusters(so.12, resolution = 0.3)

so.12$manual.cluster[which(so.12$seurat_clusters %in% c("0", "1"))] <- "CLL"
so.12$manual.cluster[which(so.12$seurat_clusters %in% c("2"))] <- "CLL"
so.12$manual.cluster[which(so.12$seurat_clusters %in% c("3"))] <- "CLL"
so.12$manual.cluster[which(so.12$seurat_clusters %in% c("4"))] <- "CD8 T cell"
so.12$manual.cluster[which(so.12$seurat_clusters %in% c("5"))] <- "CD4 T cell"
so.12$manual.cluster[which(so.12$seurat_clusters %in% c("6"))] <- "Cluster 6"
so.12$manual.cluster[which(so.12$seurat_clusters %in% c("7"))] <- "NK"
so.12$manual.cluster[which(so.12$seurat_clusters %in% c("8"))] <- "Myeloid"

so.12$sample <- "none"
so.12$sample[annotation$bc] <- annotation$sample

saveRDS("./data/coevolution/objects/20220504_CLL34_filtered.rds", object = so.12)
