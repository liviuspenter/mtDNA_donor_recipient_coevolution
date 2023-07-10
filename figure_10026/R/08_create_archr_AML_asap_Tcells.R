# create ArchR object of T cells with ASAP-seq data

library(ArchR)
library(gplots)
library(grid)
library(parallel)
library(Seurat)

TSB.so = readRDS('./data/10026/objects/20210601_TSB.so.rds')
TSB.so = RenameCells(TSB.so, new.names = paste0(colnames(TSB.so), '-1'))

AML.asap = loadArchRProject('./data/10026/AML.asap/')

AML.Tcell = subsetArchRProject(ArchRProj = AML.asap, cells = colnames(TSB.so)[which(TSB.so$seurat_clusters %in% c('7','8'))], 
                               outputDirectory = './data/10026/AML.Tcell', threads = 12, force = T)

AML.Tcell = addClusters(input = AML.Tcell, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 1, force = T)
set.seed(1987)
AML.Tcell = addUMAP(ArchRProj = AML.Tcell, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
AML.Tcell = addImputeWeights(AML.Tcell)

ADT.data = readRDS('./data/10026/objects/20210601_ADT.rds')
rownames(ADT.data) = paste0(rownames(ADT.data), '-1')
ADT.Tcell = ADT.data[AML.Tcell$cellNames,]

TSB.Tcell = TSB.so[,which(colnames(TSB.so) %in% AML.Tcell$cellNames)]
TSB.Tcell = RunPCA(TSB.Tcell)
TSB.Tcell = RunUMAP(TSB.Tcell, dims = 1:10)
TSB.Tcell = FindNeighbors(TSB.Tcell)
TSB.Tcell = FindClusters(TSB.Tcell, resolution = 1.5)
TSB.Tcell$manual.cluster = 'none'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('0'))] = 'CD8'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('1'))] = 'CD4.CM'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('2'))] = 'CD4.naive'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('3'))] = 'CD4.TEMRA'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('4','5'))] = 'CD8.senescent'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('6'))] = 'CD4.naive'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('7'))] = 'CD4.senescent'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('8'))] = 'CD4.CM'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('9'))] = 'CD4.EM'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('10'))] = 'NK'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('11'))] = 'CD4.CM'
TSB.Tcell$manual.cluster[which(TSB.Tcell$seurat_clusters %in% c('12'))] = 'CD8.CM'
TSB.Tcell$manual.cluster = factor(TSB.Tcell$manual.cluster, levels = names(Tcell.colors))
saveRDS(TSB.Tcell, './data/10026/objects/20210603_TSB_Tcell.rds')
TSB.Tcell = readRDS('./data/10026/objects/20210603_TSB_Tcell.rds')

source('./figure_10026/R/Tcell.colors.R')

DoHeatmap(TSB.Tcell, disp.min = -2, disp.max = 2, 
          features = c('CD3', 'CD4', 'CD8', 'CD62L', 'CCR7', 'CD45RO', 'CD45RA', 'CD127','CD25', 'CD28', 'CD57', 'PD1','CD39','CD38','CD56')) + 
  scale_fill_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius', type='continuous')) 

p=DoHeatmap(TSB.Tcell, disp.min = -2, disp.max = 2, group.by = 'manual.cluster', group.colors = Tcell.colors,label = F,
          features = c('CD3', 'CD4', 'CD8', 'CD62L', 'CCR7', 'CD45RO', 'CD45RA', 'CD127','CD25', 'CD28', 'CD57', 'PD1','CD39','CD38','CD56')) + 
  scale_fill_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius', type='continuous')) +
  scale_color_manual(values = Tcell.colors) + NoLegend() 
ggsave('./figure_10026/figures/heatmaps/20210603_Tcell_phenotypes.svg', width = 3, height = 4, plot = p)

svglite::svglite('./figure_10026/figures/heatmaps/brewer_celsius.svg', width = 2, height = 0.5)
BuenColors::jdb_palette(name = 'brewer_celsius', type='continuous')
dev.off()

p=DimPlot(TSB.Tcell, group.by = 'manual.cluster', cols = Tcell.colors, pt.size = 0.5) + theme(plot.title = element_blank()) + NoLegend() + NoAxes()
ggsave('./figure_10026/figures/umaps/20210603_Tcells_ASAP_seq.png', width = 4, height = 4, dpi = 600)

AML.Tcell = addCellColData(AML.Tcell, data = as.character(TSB.Tcell$seurat_clusters), name = 'seurat_clusters', cells = colnames(TSB.Tcell))
AML.Tcell = addCellColData(AML.Tcell, data = as.character(TSB.Tcell$manual.cluster), name = 'manual_clusters', cells = colnames(TSB.Tcell))

AML.Tcell <- addGroupCoverages(AML.Tcell, groupBy = "manual_clusters", threads = 12)
AML.Tcell = addReproduciblePeakSet(AML.Tcell, groupBy = 'manual_clusters', threads = 12)
AML.Tcell = addPeakMatrix(AML.Tcell, threads = 12)
AML.Tcell <- addMotifAnnotations(AML.Tcell, motifSet = "cisbp", name = "Motif")
AML.Tcell <- addBgdPeaks(AML.Tcell)
AML.Tcell <- addDeviationsMatrix(AML.Tcell, peakAnnotation = "Motif", threads = 12, force = T)
saveArchRProject(AML.Tcell)
AML.Tcell = loadArchRProject('./data/10026/AML.Tcell/')

p=plotEmbedding(AML.Tcell, name = 'manual_clusters', pal = Tcell.colors)
plotPDF(plotList = list(p), ArchRProj = AML.Tcell, name = 'Tcell.clusters', width = 4, height = 4, addDOC = F)

p = plotBrowserTrack(AML.Tcell, groupBy = 'manual_clusters', geneSymbol = c('TCF7'), upstream = 40000, downstream = 40000,
                     pal = Tcell.colors, useGroups = names(Tcell.colors))
svglite::svglite('./figure_10026/figures/browsertracks/20210603_Tcells_TCF7.svg', width = 4, height = 4)
grid.draw(p$TCF7)
dev.off()

p = plotBrowserTrack(AML.Tcell, groupBy = 'manual_clusters', geneSymbol = c('TOX'), upstream = 450000, downstream = 100000,
                     pal = Tcell.colors, useGroups = names(Tcell.colors))
svglite::svglite('./figure_10026/figures/browsertracks/20210603_Tcells_TOX.svg', width = 4, height = 4)
grid.draw(p$TOX)
dev.off()
