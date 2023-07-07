library(ArchR)
library(gplots)
library(Seurat)
library(ComplexHeatmap)

MART1.mito = loadArchRProject('/Users/shaka87/ArchRProjects/MART1/MART1.mito/')

setwd('/Users/liviuspenter/dfci/asap_seq/')

row.order = c('Tcrb','CD3', 'CD4', 'CD8', 'CD62L', 'CCR7', 'CD45RO', 'CD45RA', 'CD127','CD25', 'CD28', 'CD57', 'PD1','CD39','CD38','CD56',
              'CD14', 'CD16','CD33', 'CD11c', 'CD117', 'CD19', 'CD138')

# subset TSB data to cells in MART1.mito object and rerun UMAP + clustering
TSB.so = readRDS('./data/objects/20211031_MART1_TSB.so.rds')
TSB.so = TSB.so[,which(colnames(TSB.so) %in% MART1.mito$cellNames)]
TSB.so = RunUMAP(TSB.so, dims = 1:20)
TSB.so = FindNeighbors(TSB.so)
TSB.so = FindClusters(TSB.so, resolution = 1)

DoHeatmap(TSB.so, disp.min = -2, disp.max = 2, features = row.order) + 
  scale_fill_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius', type='continuous')) 

TSB.so$manual.cluster = 'none'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('0'))] = 'Mono'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('1'))] = 'CD4.memory'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('2'))] = 'CD8.naive'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('3'))] = 'CD4.naive'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('4'))] = 'CD8.EM'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('5'))] = 'CD8.TEMRA'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('6'))] = 'CD4.CD57hi'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('7'))] = 'unspecific'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('8'))] = 'B cell'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('9'))] = 'Mono'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('10'))] = 'CD4.activated'
TSB.so$manual.cluster[which(TSB.so$seurat_clusters %in% c('11'))] = 'NK'
saveRDS(TSB.so, file = './data/objects/20211031_MART1_TSB_mito.so.rds')

# add donor status to TSB object
boo = data.frame(cellName = MART1.mito$cellNames, 
                 donor = MART1.mito$donor)
rownames(boo) = boo$cellName
boo = boo[colnames(TSB.so),]
TSB.so$donor = boo$donor
saveRDS(TSB.so, file = './data/objects/20211031_MART1_TSB_mito.so.rds')

df = as.data.frame(t(GetAssayData(TSB.so)[c('CD3', 'CD38','CD127'),]))
df$donor = TSB.so$donor
ggplot() + 
  geom_density_2d(data=df, aes(x=CD3, y=CD38)) +
  geom_point(data=df[which(df$donor == '108'),], aes(x=CD3, y=CD38), size=0.5) +
  theme_classic()


MART1.cells.color = c('Mono' = 'darkgreen', 
                      'B cell' = 'red', 
                      'CD8.naive' = RColorBrewer::brewer.pal(name = 'Purples', n=4)[2],
                      'CD8.EM' = RColorBrewer::brewer.pal(name = 'Purples', n=4)[3],
                      'CD8.TEMRA' = RColorBrewer::brewer.pal(name = 'Purples', n=4)[4],
                      'CD4.naive' = RColorBrewer::brewer.pal(name = 'Blues', n=5)[2],
                      'CD4.memory' = RColorBrewer::brewer.pal(name = 'Blues', n=5)[3],
                      'CD4.CD57hi' = RColorBrewer::brewer.pal(name = 'Blues', n=5)[4],
                      'CD4.activated' = RColorBrewer::brewer.pal(name = 'Blues', n=5)[5],
                      'NK' = 'black',
                      'Tcrb.pos.PD1.hi' = 'grey',
                      'none' = 'grey')
TSB.so$manual.cluster = factor(TSB.so$manual.cluster, levels = names(MART1.cells.color))

boo = subset(TSB.so, donor %in% c('5', '108'))
df = GetAssayData(boo)
df = t(scale(t(df)))
ha = columnAnnotation(Donor = boo$donor,
                      #manual.cluster = TSB.so$manual.cluster[which(TSB.so$donor %in% c('5', '108'))], 
                      col = list(Donor = c('5' = 'blue', '108' = 'orange')), border = T)
col_fun = circlize::colorRamp2(breaks = seq(-2,2,4/8), colors = BuenColors::jdb_palette(name='brewer_yes', n=9))
svglite::svglite('./figures/heatmaps/20211031_MART1_donor_recipient.svg', width = 6, height = 4)
Heatmap(df[row.order,], top_annotation = ha, use_raster = T, raster_quality = 10, 
        cluster_rows = F, cluster_columns = F, show_column_names = F, col = col_fun, 
        column_title_rot = 90, column_title_side = 'bottom', column_title_gp = gpar(fontsize=10), 
        row_title_gp = gpar(fontsize=10), row_names_side = 'left',
        column_split = boo$manual.cluster, border = T)
dev.off()

MART1.mito = addCellColData(MART1.mito, data = as.character(TSB.so$manual.cluster), name = 'manual.cluster', cells = colnames(TSB.so), force=T)
saveArchRProject(MART1.mito)
p=plotEmbedding(MART1.mito, name = 'donor', pal = c('108' = 'orange', '5' = 'purple'))
plotPDF(list(p), name='MART1.donor', ArchRProj = MART1.mito, addDOC = F)

p=plotEmbedding(MART1.mito, name = 'manual.cluster', pal=MART1.cells.color)
plotPDF(list(p), name='MART1.celltypes', ArchRProj = MART1.mito, addDOC = F)

# difference between both donor-derived clusters

MART1.mito$boo = 'none'
MART1.mito$boo[which(MART1.mito$Clusters == 'C3' & MART1.mito$donor == '108')] = '108.1'
MART1.mito$boo[which(MART1.mito$Clusters == 'C7' & MART1.mito$donor == '108')] = '108.2'

markersGS <- getMarkerFeatures(
  ArchRProj = MART1.mito, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "boo",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = c('108.1'), 
  bgdGroups = c('108.2')
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")

