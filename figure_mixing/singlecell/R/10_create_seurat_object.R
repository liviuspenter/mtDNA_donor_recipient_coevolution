setwd('/Users/shaka87/dfci/asap_seq/')

library(ArchR)
library(ggplot2)
library(Seurat)

seurat.data = Read10X(data.dir = paste0('./data/artifical_mixing/Pool91-1_22/filtered_feature_bc_matrix/'))
so1 <- CreateSeuratObject(counts = seurat.data, project = 'CLL1', min.cells = 3, min.features = 200)
so1[['percent.mt']] <- PercentageFeatureSet(so1, pattern = '^MT-')
so1 = RenameCells(so1, add.cell.id = 'CLL1')

seurat.data = Read10X(data.dir = paste0('./data/artifical_mixing/Pool91-1_24/filtered_feature_bc_matrix/'))
so3 <- CreateSeuratObject(counts = seurat.data, project = 'CLL3', min.cells = 3, min.features = 200)
so3[['percent.mt']] <- PercentageFeatureSet(so3, pattern = '^MT-')
so3 = RenameCells(so3, add.cell.id = 'CLL3')

CLL.mix.so = merge(so1, list(so3) )

CLL.mix.so <- subset(CLL.mix.so, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 20)

CLL.mix.so = NormalizeData(CLL.mix.so)
CLL.mix.so = FindVariableFeatures(CLL.mix.so)
CLL.mix.so = ScaleData(CLL.mix.so)
CLL.mix.so = RunPCA(CLL.mix.so)
CLL.mix.so = RunUMAP(CLL.mix.so, dims = 1:30)
CLL.mix.so = FindNeighbors(CLL.mix.so)
CLL.mix.so = FindClusters(CLL.mix.so, resolution = 0.3)

CLL.mix.so$manual.cluster = 'CLL1'
CLL.mix.so$manual.cluster[which(CLL.mix.so$seurat_clusters %in% c('3', '5'))] = 'T cell'
CLL.mix.so$manual.cluster[which(CLL.mix.so$seurat_clusters %in% c('4'))] = 'Mono'
CLL.mix.so$manual.cluster[which(CLL.mix.so$seurat_clusters %in% c('0'))] = 'CLL3'

saveRDS('./data/artifical_mixing/CLL_mix_so.rds', object = CLL.mix.so)