library(ArchR)
library(gplots)
library(grid)
library(Seurat)

addArchRGenome('hg38')

setwd('/Users/shaka87/dfci/asap_seq/')

TSB.so = readRDS('./data/objects/20210601_TSB.so.rds')
AML.1012.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1012.mito/')
AML.1012.mito$manual.clusters = c('none')
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% 'C9')] = c('CD4')
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% 'C10')] = c('CD8')
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c('C7'))] = c('erythroid')
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c('C4'))] = c('HSCT')
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c('C5','C6'))] = c('GMP')
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c('C2'))] = c('Mono1')
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c('C3'))] = c('Mono2')
AML.1012.mito$manual.clusters[which(AML.1012.mito$Clusters %in% c('C1'))] = c('DC')

AML.1012.mito.colors = c('HSCT' = RColorBrewer::brewer.pal(n=6, name='Reds')[6], 
                         'GMP' = RColorBrewer::brewer.pal(n=6, name='Reds')[5], 
                         'Mono1' = RColorBrewer::brewer.pal(n=6, name='Reds')[4],
                         'Mono2' = RColorBrewer::brewer.pal(n=6, name='Reds')[3],
                         'erythroid' = RColorBrewer::brewer.pal(n=6, name='Reds')[2], 
                         'DC' = 'yellow',
                         'CD4' = RColorBrewer::brewer.pal(n=3, name='Blues')[2],
                         'CD8' = RColorBrewer::brewer.pal(n=3, name='Blues')[3],
                         'none' = 'grey')

p=plotEmbedding(AML.1012.mito, name = 'manual.clusters', pal = c(AML.1012.mito.colors))
plotPDF(ArchRProj = AML.1012.mito, name = 'manual.clusters', plotList = list(p), addDOC = F)

combined.frequencies = readRDS('./data/mtDNA/20210429_AML1012_combined_mutation_frequencies.rds')

#combined.frequencies = combined.frequencies[-which(rownames(combined.frequencies) %in% donor.variants),]

# select mutations than mark between 500 and 1000 cells
cells.marked = apply(combined.frequencies, 1, function(x) length(which(x != 0)))
combined.frequencies = combined.frequencies[cells.marked %in% c(seq(50,1000)),]

combined.frequencies = combined.frequencies[rowMeans(as.matrix(combined.frequencies)) > 0.0025,]

cells = AML.1012.mito$cellNames
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

cells.ordered = c()
for (cluster in names(AML.1012.mito.colors)) {
  boo = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == cluster)]
  if (length(boo) > 500) {
    boo = boo[sample(length(boo), 500)]
  }
  names(boo) = rep(cluster, length(boo))
  cells.ordered = c(cells.ordered, boo)
}

# plot heatmap
ha.top = ComplexHeatmap::columnAnnotation(cluster = names(cells.ordered),
                                          show_legend=F,show_annotation_name = T, 
                                          annotation_label = c('cluster' = 'phenotype'), 
                                          col = list(cluster = AML.1012.mito.colors),
                                          annotation_name_gp = gpar(fontsize= 8) )

col_fun = circlize::colorRamp2(seq(0.0001, 0.01, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))

svglite::svglite('./figures/heatmaps/20210611_AML1012_myeloid_mtDNA.svg', width = 4, height = 4.5)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[,cells.ordered]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = T,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = paste0(nrow(combined.frequencies), ' mtDNA mutations'), row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = paste0(length(cells.ordered), ' cells'), 
                        column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, col = col_fun, 
                        row_names_gp = gpar(fontsize = 8), use_raster = F)
dev.off()

svglite::svglite('./figures/heatmaps/20210611_AML1012_myeloid_mtDNA2.svg', width = 4, height = 4.5)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[,cells.ordered]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = T, cluster_rows = T,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = paste0(nrow(combined.frequencies), ' mtDNA mutations'), row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = paste0(length(cells.ordered), ' cells'), 
                        column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, col = col_fun, 
                        row_names_gp = gpar(fontsize = 8), use_raster = F)
dev.off()

ADT.mat = GetAssayData(TSB.so, slot = 'scale.data')[,gsub(cells.ordered, pattern = '-1', replacement = '')]

col_fun = circlize::colorRamp2(seq(-2, 2, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'brewer_celsius', n=9))(100))
svglite::svglite('./figures/heatmaps/20210611_AML1012_myeloid_ADT.svg', width = 4, height = 2)
ComplexHeatmap::Heatmap(as.matrix(ADT.mat[c('CD38', 'CD33', 'CD117', 'CD14', 'CD39', 'CD11c', 'CD3', 'CD4', 'CD8', 'CD57'),]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = T,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'ASAP-seq', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = paste0(length(cells.ordered), ' cells'),  
                        column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, col = col_fun, 
                        row_names_gp = gpar(fontsize = 8), use_raster = F)
dev.off()


# cluster over time
cells.mat = data.frame(sample = as.character(), cluster = as.character(), cells = as.numeric())
for (s in c('AML1012_1', 'AML1012_2', 'AML1012_4')) {
  for (cluster in names(AML.1012.mito.colors)) {
    boo = length(which(AML.1012.mito$Sample.hash == s & AML.1012.mito$manual.clusters == cluster))
    cells.mat = rbind(cells.mat, data.frame(sample = s, cluster = cluster,                                             
                                            cells = boo,
                                            freq = boo / length(which(AML.1012.mito$Sample.hash == s))))
    
  }
}
cells.mat$cluster = factor(cells.mat$cluster, levels = names(AML.1012.mito.colors))
ggplot(data=cells.mat, aes(x=sample, y=100*freq, group=cluster, fill=cluster)) + 
  geom_col() + 
  scale_fill_manual(values = AML.1012.mito.colors) +
  scale_x_discrete(labels = c('Screening', 'EOLN', 'EOT')) + 
  scale_y_continuous('%cells') + 
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.title.x = element_blank(),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text('Arial', size=10, color='black', angle=90, hjust=1, vjust=0.5),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figures/plots/20210611_AML1012_cluster_kinetics.svg', width = 1.5, height = 2.5)

mtDNA.kinetics = names(sort(HSCT.df[3,order(colMeans(HSCT.df))], decreasing = T))

### mtDNA mutations within HSCT cluster
HSCT.1 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'HSCT' & AML.1012.mito$Sample.hash == 'AML1012_1')]
HSCT.2 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'HSCT' & AML.1012.mito$Sample.hash == 'AML1012_2')]
HSCT.4 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'HSCT' & AML.1012.mito$Sample.hash == 'AML1012_4')]

HSCT.df = data.frame(rbind(rowMeans(combined.frequencies[HSCT.1]), rowMeans(combined.frequencies[HSCT.2])))
HSCT.df = rbind(HSCT.df, rowMeans(combined.frequencies[HSCT.4]))
colnames(HSCT.df) = rownames(combined.frequencies)
rownames(HSCT.df) = c('Screening', 'EOLN', 'EOT')

col_fun = circlize::colorRamp2(seq(0.0001, 0.01, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))
svglite::svglite('./figures/heatmaps/20210612_AML1012_HSCT_mtDNA_kinetics.svg', width = 2, height = 3)
ComplexHeatmap::Heatmap(as.matrix(t(HSCT.df[,mtDNA.kinetics])),
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = T, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'HSCT mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = NA, column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        col = col_fun, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),use_raster = F)
dev.off()

### mtDNA mutations within GMP cluster
GMP.1 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'GMP' & AML.1012.mito$Sample.hash == 'AML1012_1')]
GMP.2 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'GMP' & AML.1012.mito$Sample.hash == 'AML1012_2')]
GMP.4 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'GMP' & AML.1012.mito$Sample.hash == 'AML1012_4')]

GMP.df = data.frame(rbind(rowMeans(combined.frequencies[GMP.1]), rowMeans(combined.frequencies[GMP.2])))
GMP.df = rbind(GMP.df, rowMeans(combined.frequencies[GMP.4]))
colnames(GMP.df) = rownames(combined.frequencies)
rownames(GMP.df) = c('Screening', 'EOLN', 'EOT')

col_fun = circlize::colorRamp2(seq(0.0001, 0.01, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))
svglite::svglite('./figures/heatmaps/20210612_AML1012_GMP_mtDNA_kinetics.svg', width = 2, height = 3)
ComplexHeatmap::Heatmap(as.matrix(t(GMP.df[,mtDNA.kinetics])),
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = T, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'GMP mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = NA, column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        col = col_fun, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),use_raster = F)
dev.off()

### mtDNA mutations within Mono1 cluster
Mono1.1 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'Mono1' & AML.1012.mito$Sample.hash == 'AML1012_1')]
Mono1.2 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'Mono1' & AML.1012.mito$Sample.hash == 'AML1012_2')]
Mono1.4 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'Mono1' & AML.1012.mito$Sample.hash == 'AML1012_4')]

Mono1.df = data.frame(rbind(rowMeans(combined.frequencies[Mono1.1]), rowMeans(combined.frequencies[Mono1.2])))
Mono1.df = rbind(Mono1.df, rowMeans(combined.frequencies[Mono1.4]))
colnames(Mono1.df) = rownames(combined.frequencies)
rownames(Mono1.df) = c('Screening', 'EOLN', 'EOT')

col_fun = circlize::colorRamp2(seq(0.0001, 0.01, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))
svglite::svglite('./figures/heatmaps/20210612_AML1012_Mono1_mtDNA_kinetics.svg', width = 2, height = 3)
ComplexHeatmap::Heatmap(as.matrix(t(Mono1.df[,mtDNA.kinetics])),
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = T, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'Mono1 mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = NA, column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        col = col_fun, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),use_raster = F)
dev.off()

### mtDNA mutations within Mono2 cluster
Mono2.1 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'Mono2' & AML.1012.mito$Sample.hash == 'AML1012_1')]
Mono2.2 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'Mono2' & AML.1012.mito$Sample.hash == 'AML1012_2')]
Mono2.4 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'Mono2' & AML.1012.mito$Sample.hash == 'AML1012_4')]

Mono2.df = data.frame(rbind(rowMeans(combined.frequencies[Mono2.1]), rowMeans(combined.frequencies[Mono2.2])))
Mono2.df = rbind(Mono2.df, rowMeans(combined.frequencies[Mono2.4]))
colnames(Mono2.df) = rownames(combined.frequencies)
rownames(Mono2.df) = c('Screening', 'EOLN', 'EOT')

col_fun = circlize::colorRamp2(seq(0.0001, 0.01, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))
svglite::svglite('./figures/heatmaps/20210612_AML1012_Mono2_mtDNA_kinetics.svg', width = 2, height = 3)
ComplexHeatmap::Heatmap(as.matrix(t(Mono2.df[,mtDNA.kinetics])),
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = T, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'Mono2 mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = NA, column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        col = col_fun, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),use_raster = F)
dev.off()

### mtDNA mutations within Ery cluster
Ery.1 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'erythroid' & AML.1012.mito$Sample.hash == 'AML1012_1')]
Ery.2 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'erythroid' & AML.1012.mito$Sample.hash == 'AML1012_2')]
Ery.4 = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters == 'erythroid' & AML.1012.mito$Sample.hash == 'AML1012_4')]

Ery.df = data.frame(rbind(rowMeans(combined.frequencies[Ery.1]), rowMeans(combined.frequencies[Ery.2])))
Ery.df = rbind(Ery.df, rowMeans(combined.frequencies[Ery.4]))
colnames(Ery.df) = rownames(combined.frequencies)
rownames(Ery.df) = c('Screening', 'EOLN', 'EOT')

col_fun = circlize::colorRamp2(seq(0.0001, 0.01, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))
svglite::svglite('./figures/heatmaps/20210612_AML1012_Ery_mtDNA_kinetics.svg', width = 2, height = 3)
ComplexHeatmap::Heatmap(as.matrix(t(Ery.df[,mtDNA.kinetics])),
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = T, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'Ery mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = NA, column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        col = col_fun, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),use_raster = F)
dev.off()
