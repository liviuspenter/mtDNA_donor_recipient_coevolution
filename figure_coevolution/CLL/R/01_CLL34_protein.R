# create Seurat object based on Total-seq D data

library(dplyr)
library(ggplot2)
library(Seurat)

protein.1 = as.data.frame(data.table::fread('./data/coevolution/CLL3/protein_reads.csv'))
protein.1$Barcode = paste0('CLL3#', protein.1$Barcode)
rownames(protein.1) = protein.1$Barcode
protein.1 = protein.1[,-c(1,2)]
protein.1 = t(protein.1)
protein.2 = as.data.frame(data.table::fread('./data/coevolution/CLL4/protein_reads.csv'))
protein.2$Barcode = paste0('CLL4#', protein.2$Barcode)
rownames(protein.2) = protein.2$Barcode
protein.2 = protein.2[,-c(1,2)]
protein.2 = t(protein.2)

so.1 = CreateSeuratObject(counts = protein.1, assay = 'ADT', project = 'CLL3')
so.2 = CreateSeuratObject(counts = protein.2, assay = 'ADT', project = 'CLL4')
so.12 = merge(so.1, so.2)

so.12 = NormalizeData(so.12, normalization.method = 'CLR', margin = 2)
so.12 = FindVariableFeatures(so.12)
so.12 = ScaleData(so.12, features = rownames(so.12))
so.12 = RunPCA(so.12)
so.12 = RunUMAP(so.12, dims = 1:15)
so.12 = FindNeighbors(so.12)
so.12 = FindClusters(so.12, resolution = 0.2)

pbmc.markers <- FindAllMarkers(so.12, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

DoHeatmap(subset(so.12, downsample = 200), features = c('CD5', pbmc.markers$gene)) + 
  scale_fill_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_yes'))

so.12$manual.cluster = 'CLL'
so.12$manual.cluster[which(so.12$seurat_clusters == '3')] = 'CD8 T cell'
so.12$manual.cluster[which(so.12$seurat_clusters == '4')] = 'CD4 T cell'
so.12$manual.cluster[which(so.12$seurat_clusters == '5')] = 'Myeloid'
so.12$manual.cluster[which(so.12$seurat_clusters == '6')] = 'Cluster 6'

saveRDS('./data/coevolution/objects/20220504_CLL34.rds', object = so.12)
