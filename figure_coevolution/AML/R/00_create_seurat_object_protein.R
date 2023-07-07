# read in protein expression data and create Seurat object for further processing

library(ggplot2)
library(Seurat)

# read data
protein.1 = as.data.frame(data.table::fread('./data/coevolution/1026_1/protein_reads.csv'))
protein.1$Barcode = paste0('1026.1#', protein.1$Barcode)
rownames(protein.1) = protein.1$Barcode
protein.1 = protein.1[,-c(1,2)]
protein.1 = t(protein.1)
protein.3 = as.data.frame(data.table::fread('./data/coevolution/1026_3/protein_reads.csv'))
protein.3$Barcode = paste0('1026.3#', protein.3$Barcode)
rownames(protein.3) = protein.3$Barcode
protein.3 = protein.3[,-c(1,2)]
protein.3 = t(protein.3)

# create Seurat objects and merge
so.1 = CreateSeuratObject(counts = protein.1, assay = 'ADT', project = '1026.1')
so.3 = CreateSeuratObject(counts = protein.3, assay = 'ADT', project = '1026.3')
so.13 = merge(so.1, so.3)

# standard Seurat procedures
so.13 = NormalizeData(so.13, normalization.method = 'CLR', margin = 2)
so.13 = FindVariableFeatures(so.13)
so.13 = ScaleData(so.13)
so.13 = RunPCA(so.13)
so.13 = RunUMAP(so.13, dims = 1:10)
so.13 = FindNeighbors(so.13)
so.13 = FindClusters(so.13, resolution = 0.3)

# identify markers for cluster annotation
pbmc.markers <- FindAllMarkers(so.13, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC))

# cluster annotation
so.13$manual.cluster = 'none'
so.13$manual.cluster[which(so.13$seurat_clusters == '0')] = 'Progenitor'
so.13$manual.cluster[which(so.13$seurat_clusters == '1')] = 'Mono'
so.13$manual.cluster[which(so.13$seurat_clusters == '2')] = 'HSC'
so.13$manual.cluster[which(so.13$seurat_clusters == '3')] = 'Erythroid'
so.13$manual.cluster[which(so.13$seurat_clusters == '4')] = 'CD4'
so.13$manual.cluster[which(so.13$seurat_clusters == '5')] = 'Plasma'
so.13$manual.cluster[which(so.13$seurat_clusters == '6')] = 'CD8'
so.13$manual.cluster = factor(so.13$manual.cluster, levels = c('HSC', 'Progenitor', 'Mono', 'Erythroid', 'CD4', 'CD8', 'Plasma'))

protein.order = c('CD34', 'CD117','CD123', 'CD38', 'HLA-DR', 'CD56','CD33', 'CD64','CD14', 'CD16', 'CD10', 'CD11b', 'CD11c', 'CD45', 'CD71', 'CD141','CD49d',
                  'CD2', 'CD5', 'CD7','CD3', 'CD4', 'CD8', 'CD45RA', 'CD45RO', 'CD62L', 'CD25','CD138', 'CD19', 'CD22', 'CD30')

cluster.colors = c('HSC' = RColorBrewer::brewer.pal(name = 'YlOrRd', n = 5)[5],
                   'Progenitor' = RColorBrewer::brewer.pal(name = 'YlOrRd', n = 5)[4],
                   'Mono' = RColorBrewer::brewer.pal(name = 'YlOrRd', n = 5)[3],
                   'Erythroid' = RColorBrewer::brewer.pal(name = 'YlOrRd', n = 5)[2],
                   'CD4' = 'lightblue', 
                   'CD8' = 'darkblue', 
                   'Plasma' = 'darkgreen')

# heatmap of lineage markers
p=DoHeatmap(subset(so.13, downsample=200), group.by = 'manual.cluster', features = protein.order, group.colors = cluster.colors, label = F) + 
  scale_fill_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_yes'), na.value = 'white') +
  NoLegend() +
  theme(axis.text.y = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/AML/figures/20220112_AML1026_cluster_features.png', width = 2.5, height = 3.5, plot = p, dpi = 600)

# UMAP of major celltypes
p=DimPlot(so.13, group.by = 'manual.cluster', cols = cluster.colors) + NoLegend() + NoAxes() + 
  theme(plot.title = element_blank())
ggsave('./figure_coevolution/AML/figures/20220112_AML1026_UMAP.png', width = 3, height = 3, dpi = 600, plot = p)

# UMAP of samples
p=DimPlot(so.13, group.by = 'orig.ident', cols = c('blue', 'yellow')) + NoLegend() + NoAxes() + 
  theme(plot.title = element_blank())
ggsave('./figure_coevolution/AML/figures/20220112_AML1026_UMAP_samples.png', width = 3, height = 3, dpi = 600, plot = p)

# kinetics of major celltypes
df = as.data.frame(prop.table(table(Idents(so.13), so.13$orig.ident), margin = 2))
ggplot(data=df, aes(x=Var2, y=100*Freq, fill=Var1)) + geom_col() + scale_fill_manual(values = cluster.colors) +
  scale_x_discrete(breaks = c('1026.1', '1026.3'), labels = c('Screening', 'On treatment')) + 
  scale_y_continuous('% cells') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_coevolution/AML/figures/20220112_AML1026_cluster_kinetics.svg', width = 1.2, height = 2.5)

# brewer_yes color key
svglite::svglite('./figure_coevolution/AML/figures/legend_brewer_yes.svg', width = 2, height = 0.5)
BuenColors::jdb_palette(name = 'brewer_yes', type = 'continuous')
dev.off()

# solar_rojos color key
svglite::svglite('./figure_coevolution/AML/figures/legend_solar_rojos.svg', width = 2, height = 0.5)
BuenColors::jdb_palette(name = 'solar_rojos', type = 'continuous')
dev.off()

saveRDS('./data/coevolution/objects/20220112_AML1026.rds', object = so.13)

# generate Seurat object only for T cells
so.13.Tcell = subset(so.13,manual.cluster %in% c('CD4', 'CD8'))
so.13.Tcell = RunUMAP(so.13.Tcell, dims = 1:10)
so.13.Tcell = FindNeighbors(so.13.Tcell)
so.13.Tcell = FindClusters(so.13.Tcell, resolution = 2)

# annotate
so.13.Tcell$manual.cluster = 'none'
so.13.Tcell$manual.cluster[which(so.13.Tcell$seurat_clusters == 0)] = 'NK'
so.13.Tcell$manual.cluster[which(so.13.Tcell$seurat_clusters == 1)] = 'CD8 memory'
so.13.Tcell$manual.cluster[which(so.13.Tcell$seurat_clusters == 2)] = 'CD4 CM activated'
so.13.Tcell$manual.cluster[which(so.13.Tcell$seurat_clusters == 3)] = 'CD4 CM'
so.13.Tcell$manual.cluster[which(so.13.Tcell$seurat_clusters == 4)] = 'CD4 naive'
so.13.Tcell$manual.cluster[which(so.13.Tcell$seurat_clusters == 5)] = 'CD4 EM'
so.13.Tcell$manual.cluster[which(so.13.Tcell$seurat_clusters == 6)] = 'CD8'
so.13.Tcell$manual.cluster[which(so.13.Tcell$seurat_clusters == 7)] = 'CD4 CM'
so.13.Tcell$manual.cluster[which(so.13.Tcell$seurat_clusters == 8)] = 'CD4 naive'
so.13.Tcell$manual.cluster[which(so.13.Tcell$seurat_clusters == 9)] = 'CD4'
so.13.Tcell$manual.cluster = factor(so.13.Tcell$manual.cluster, levels = c('CD4', 'CD4 naive', 'CD4 CM', 'CD4 CM activated', 'CD4 EM', 'CD8', 'CD8 memory', 'NK'))
saveRDS('./data/coevolution/objects/20220112_AML1026_Tcell.rds', object = so.13.Tcell)

# heatmap of T cell subpopulations 
cluster.colors.T = c('CD4' = RColorBrewer::brewer.pal(name = 'Blues', n = 6)[2],
                     'CD4 naive' = RColorBrewer::brewer.pal(name = 'Blues', n = 6)[3],
                     'CD4 CM' = RColorBrewer::brewer.pal(name = 'Blues', n = 6)[4],
                     'CD4 CM activated' = RColorBrewer::brewer.pal(name = 'Blues', n = 6)[5],
                     'CD4 EM' = RColorBrewer::brewer.pal(name = 'Blues', n = 6)[6], 
                     'CD8' = RColorBrewer::brewer.pal(name = 'Reds', n = 3)[2],
                     'CD8 memory' = RColorBrewer::brewer.pal(name = 'Reds', n = 3)[3],
                     'NK' = 'black')

DoHeatmap(so.13.Tcell, group.by = 'manual.cluster', features = c('CD3', 'CD56', 'CD16','CD4', 'CD8', 'CD62L', 'CD45RA', 'CD45RO', 'CD44', 'CD69'),
          group.colors = cluster.colors.T, label = F) + 
  scale_fill_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_yes'), na.value = 'white') +
  NoLegend() +
  theme(axis.text.y = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/AML/figures/20220113_AML1026_Tcell_phenotypes.png', width = 2, height = 2, dpi = 600)