library(ArchR)
library(gplots)
library(grid)
library(Seurat)

addArchRGenome('hg38')

setwd('/Users/shaka87/dfci/asap_seq/')

TSB.so = readRDS('./data/objects/20210601_TSB.so.rds')
AML.1010.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1010.mito/')
AML.1010.mito$manual.clusters = c('AML')
AML.1010.mito$manual.clusters[which(AML.1010.mito$Clusters %in% 'C3')] = c('CD4')
AML.1010.mito$manual.clusters[which(AML.1010.mito$Clusters %in% 'C2')] = c('CD8')
AML.1010.mito$manual.clusters[which(AML.1010.mito$Clusters %in% c('C1','C6'))] = c('erythroid')
AML.1010.mito$manual.clusters[which(AML.1010.mito$Clusters %in% c('C7','C8'))] = c('HSCT')
AML.1010.mito$manual.clusters[which(AML.1010.mito$Clusters %in% c('C5'))] = c('GMP')
AML.1010.mito$manual.clusters[which(AML.1010.mito$Clusters %in% c('C4'))] = c('Mono')

AML.1010.mito.colors = c('HSCT' = RColorBrewer::brewer.pal(n=5, name='Reds')[5], 
                         'GMP' = RColorBrewer::brewer.pal(n=5, name='Reds')[4], 
                         'Mono' = RColorBrewer::brewer.pal(n=5, name='Reds')[3],
                         'erythroid' = RColorBrewer::brewer.pal(n=5, name='Reds')[2], 
                         'CD4' = RColorBrewer::brewer.pal(n=3, name='Blues')[2],
                         'CD8' = RColorBrewer::brewer.pal(n=3, name='Blues')[3])

p=plotEmbedding(AML.1010.mito, name = 'manual.clusters', pal = c(AML.1010.mito.colors))
plotPDF(ArchRProj = AML.1010.mito, name = 'manual.clusters', plotList = list(p), addDOC = F)

donor.variants = gtools::mixedsort(c('195T>C', '150C>T', '6146A>G', '1811A>G', '10907T>C', '6047A>G', '5999T>C', '14866C>T', '11009T>C', '11467A>G', '9070T>G', '15693T>C', '12308A>G',
                                     '12372G>A', '4646T>C', '14620C>T', '15530T>C', '4811A>G', '499G>A', '16356T>C', '11332C>T', '16179C>T'))
recipient.variants = gtools::mixedsort(c('16126T>C', '4216T>C', '10084T>C', '489T>C', '462C>T', '11251A>G', '12612A>G', '15452C>A', '3010G>A', '14798T>C', '13708G>A', 
                                         '16069C>T', '16319G>A', '295C>T', '55T>C', '56A>G', '185G>A', '228G>A'))

combined.frequencies = readRDS('./data/mtDNA/20210429_AML1010_combined_mutation_frequencies.rds')
combined.frequencies = combined.frequencies[-which(rownames(combined.frequencies) == '310T>C'),]
combined.frequencies = combined.frequencies[-which(rownames(combined.frequencies) %in% donor.variants),]

# select mutations than mark between 500 and 5000 cells
cells.marked = apply(combined.frequencies, 1, function(x) length(which(x != 0)))
combined.frequencies = combined.frequencies[cells.marked %in% c(seq(50,5000)),]

combined.frequencies = combined.frequencies[rowVars(as.matrix(combined.frequencies)) > 0.0005,]

cells = AML.1010.mito$cellNames
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

cells.ordered = c()
for (cluster in names(AML.1010.mito.colors)) {
  boo = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == cluster)]
  if (length(boo) > 1000) {
    boo = boo[sample(length(boo), 1000)]
  }
  names(boo) = rep(cluster, length(boo))
  cells.ordered = c(cells.ordered, boo)
}

# plot heatmap
ha.top = ComplexHeatmap::columnAnnotation(cluster = names(cells.ordered),
                                          show_legend=F,show_annotation_name = T, 
                                          annotation_label = c('cluster' = 'phenotype'), 
                                          col = list(cluster = AML.1010.mito.colors),
                                          annotation_name_gp = gpar(fontsize= 8) )

col_fun = circlize::colorRamp2(seq(0.001, 0.1, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))

svglite::svglite('./figures/heatmaps/20210611_AML1010_myeloid_mtDNA.svg', width = 4, height = 4.5)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[,cells.ordered]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = T,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = '28 mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = '4237 cells', column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, col = col_fun, 
                        row_names_gp = gpar(fontsize = 8), use_raster = F)
dev.off()

svglite::svglite('./figures/heatmaps/20210611_AML1010_myeloid_mtDNA2.svg', width = 4, height = 4.5)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[,cells.ordered]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = T, cluster_rows = T,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = '28 mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = '4237 cells', column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, col = col_fun, 
                        row_names_gp = gpar(fontsize = 8), use_raster = F)
dev.off()

ADT.mat = GetAssayData(TSB.so, slot = 'scale.data')[,gsub(cells.ordered, pattern = '-1', replacement = '')]

col_fun = circlize::colorRamp2(seq(-2, 2, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'brewer_celsius', n=9))(100))
svglite::svglite('./figures/heatmaps/20210611_AML1010_myeloid_ADT.svg', width = 4, height = 2)
ComplexHeatmap::Heatmap(as.matrix(ADT.mat[c('CD38', 'CD33', 'CD117', 'CD14', 'CD39', 'CD11c', 'CD3', 'CD4', 'CD8', 'CD57'),]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = T,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'ASAP-seq', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = '4237 cells', column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, col = col_fun, 
                        row_names_gp = gpar(fontsize = 8), use_raster = F)
dev.off()


# cluster over time
cells.mat = data.frame(sample = as.character(), cluster = as.character(), cells = as.numeric())
for (s in c('AML1010_1', 'AML1010_2', 'AML1010_3', 'AML1010_4', 'AML1010_5')) {
  for (cluster in names(AML.1010.mito.colors)) {
    boo = length(which(AML.1010.mito$Sample.hash == s & AML.1010.mito$manual.clusters == cluster))
    cells.mat = rbind(cells.mat, data.frame(sample = s, cluster = cluster,                                             
                                            cells = boo,
                                            freq = boo / length(which(AML.1010.mito$Sample.hash == s))))
    
  }
}
cells.mat$cluster = factor(cells.mat$cluster, levels = names(AML.1010.mito.colors))
ggplot(data=cells.mat, aes(x=sample, y=100*freq, group=cluster, fill=cluster)) + 
  geom_col() + 
  scale_fill_manual(values = AML.1010.mito.colors) +
  scale_x_discrete(labels = c('Screening', 'EOLN', 'C1', 'C4', 'C10')) + 
  scale_y_continuous('%cells') + 
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.title.x = element_blank(),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text('Arial', size=10, color='black', angle=90, hjust=1, vjust=0.5),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figures/plots/20210611_AML1010_cluster_kinetics.svg', width = 1.5, height = 2)

svglite::svglite('./figures/heatmaps/solar_rojos.svg', width = 2, height = 0.5)
BuenColors::jdb_palette(name = 'solar_rojos', type = 'continuous')
dev.off()

#mtDNA.kinetics = c('8902G>C', '9525G>C', '2542G>C', '13210G>C', '9947G>C', '1782G>C',  '9655G>C', '9820G>C', '9222C>T', '139T>C', '1816G>A')
# filter mtDNA kinetics to show
mtDNA.kinetics = names(sort(HSCT.df[1,order(colMeans(HSCT.df), decreasing = T)[1:15]] / HSCT.df[5,order(colMeans(HSCT.df), decreasing = T)[1:15]]))
mtDNA.kinetics = names(sort(HSCT.df[5, mtDNA.kinetics], decreasing = T))

### mtDNA mutations within HSCT cluster
HSCT.1 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'HSCT' & AML.1010.mito$Sample.hash == 'AML1010_1')]
HSCT.2 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'HSCT' & AML.1010.mito$Sample.hash == 'AML1010_2')]
HSCT.3 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'HSCT' & AML.1010.mito$Sample.hash == 'AML1010_3')]
HSCT.4 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'HSCT' & AML.1010.mito$Sample.hash == 'AML1010_4')]
HSCT.5 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'HSCT' & AML.1010.mito$Sample.hash == 'AML1010_5')]

HSCT.df = data.frame(rbind(rowMeans(combined.frequencies[HSCT.1]), rowMeans(combined.frequencies[HSCT.2])))
HSCT.df = rbind(HSCT.df, rowMeans(combined.frequencies[HSCT.3]))
HSCT.df = rbind(HSCT.df, rowMeans(combined.frequencies[HSCT.4]))
HSCT.df = rbind(HSCT.df, rowMeans(combined.frequencies[HSCT.5]))
colnames(HSCT.df) = rownames(combined.frequencies)
rownames(HSCT.df) = c('Screening', 'EOLN', 'C1', 'C4', 'C10')

col_fun = circlize::colorRamp2(seq(0.001, 0.1, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))
svglite::svglite('./figures/heatmaps/20210612_AML1010_HSCT_mtDNA_kinetics.svg', width = 2, height = 2.5)
ComplexHeatmap::Heatmap(as.matrix(t(HSCT.df[,mtDNA.kinetics])),
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = T, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'HSCT mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = NA, column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        col = col_fun, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),use_raster = F)
dev.off()

### mtDNA mutations within GMP cluster
GMP.1 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'GMP' & AML.1010.mito$Sample.hash == 'AML1010_1')]
GMP.2 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'GMP' & AML.1010.mito$Sample.hash == 'AML1010_2')]
GMP.3 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'GMP' & AML.1010.mito$Sample.hash == 'AML1010_3')]
GMP.4 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'GMP' & AML.1010.mito$Sample.hash == 'AML1010_4')]
GMP.5 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'GMP' & AML.1010.mito$Sample.hash == 'AML1010_5')]

GMP.df = data.frame(rbind(rowMeans(combined.frequencies[GMP.1]), rowMeans(combined.frequencies[GMP.2])))
GMP.df = rbind(GMP.df, rowMeans(combined.frequencies[GMP.3]))
GMP.df = rbind(GMP.df, rowMeans(combined.frequencies[GMP.4]))
GMP.df = rbind(GMP.df, rowMeans(combined.frequencies[GMP.5]))
colnames(GMP.df) = rownames(combined.frequencies)
rownames(GMP.df) = c('Screening', 'EOLN', 'C1', 'C4', 'C10')

col_fun = circlize::colorRamp2(seq(0.001, 0.1, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))
svglite::svglite('./figures/heatmaps/20210612_AML1010_GMP_mtDNA_kinetics.svg', width = 2, height = 2.5)
ComplexHeatmap::Heatmap(as.matrix(t(GMP.df[,mtDNA.kinetics])),
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = T, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'GMP mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = NA, column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        col = col_fun, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), use_raster = F)
dev.off()

### mtDNA mutations within Ery cluster
Ery.1 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'erythroid' & AML.1010.mito$Sample.hash == 'AML1010_1')]
Ery.2 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'erythroid' & AML.1010.mito$Sample.hash == 'AML1010_2')]
Ery.3 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'erythroid' & AML.1010.mito$Sample.hash == 'AML1010_3')]
Ery.4 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'erythroid' & AML.1010.mito$Sample.hash == 'AML1010_4')]
Ery.5 = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == 'erythroid' & AML.1010.mito$Sample.hash == 'AML1010_5')]

Ery.df = data.frame(rbind(rowMeans(combined.frequencies[Ery.1]), rowMeans(combined.frequencies[Ery.2])))
Ery.df = rbind(Ery.df, rowMeans(combined.frequencies[Ery.3]))
Ery.df = rbind(Ery.df, rowMeans(combined.frequencies[Ery.4]))
Ery.df = rbind(Ery.df, rowMeans(combined.frequencies[Ery.5]))
colnames(Ery.df) = rownames(combined.frequencies)
rownames(Ery.df) = c('Screening', 'EOLN', 'C1', 'C4', 'C10')

col_fun = circlize::colorRamp2(seq(0.001, 0.1, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))
svglite::svglite('./figures/heatmaps/20210612_AML1010_Ery_mtDNA_kinetics.svg', width = 2, height = 2.5)
ComplexHeatmap::Heatmap(as.matrix(t(Ery.df[,mtDNA.kinetics])),
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = T, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'Ery mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = NA, column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        col = col_fun, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), use_raster = F)
dev.off()

### cluster by recipient and donor variants

combined.frequencies = readRDS('./data/mtDNA/20210429_AML1010_combined_mutation_frequencies.rds')
combined.frequencies = combined.frequencies[-which(rownames(combined.frequencies) == '310T>C'),]
combined.frequencies = combined.frequencies[c(recipient.variants, donor.variants),]

cells = AML.1010.mito$cellNames
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

cells.ordered = c()
for (cluster in names(AML.1010.mito.colors)) {
  boo = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == cluster)]
  if (length(boo) > 1000) {
    boo = boo[sample(length(boo), 1000)]
  }
  names(boo) = rep(cluster, length(boo))
  cells.ordered = c(cells.ordered, boo)
}

# plot heatmap
ha.top = ComplexHeatmap::columnAnnotation(cluster = names(cells.ordered),
                                          show_legend=F,show_annotation_name = T, 
                                          annotation_label = c('cluster' = 'phenotype'), 
                                          col = list(cluster = AML.1010.mito.colors),
                                          annotation_name_gp = gpar(fontsize= 8) )

col_fun = circlize::colorRamp2(seq(0.001, 0.1, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))

svglite::svglite('./figures/heatmaps/20210624_AML1010_myeloid_mtDNA_recipient_donor_variants.svg', width = 4, height = 4.5)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[,cells.ordered]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = '28 mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = '4237 cells', column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, col = col_fun, 
                        row_names_gp = gpar(fontsize = 8), use_raster = F)
dev.off()

# score calculated from donor and recipient variants
chimerism.score = as.numeric(colMeans(combined.frequencies[donor.variants,cells.ordered])) - 
  as.numeric(colMeans(combined.frequencies[recipient.variants,cells.ordered]))

chimerism.score = as.numeric(colMeans(combined.frequencies[donor.variants,])) - 
  as.numeric(colMeans(combined.frequencies[recipient.variants,]))

chimerism.analysis = ifelse(chimerism.score > 0, 'donor', 'recipient')

donor.cells = cells.ordered[which(chimerism.analysis == 'donor' & !names(cells.ordered) %in% c('CD4', 'CD8'))]
recipient.cells = cells.ordered[which(chimerism.analysis == 'recipient'& !names(cells.ordered) %in% c('CD4', 'CD8'))]
