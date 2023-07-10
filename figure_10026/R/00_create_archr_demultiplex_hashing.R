# create ArchR object with all AML samples
# perform deconvolution using hashtags

library(ArchR)
library(Rsamtools)
library(gplots)
library(Seurat)

# fix for ArchR 
source('./R/20200708_filter_doublets_archr_fix.R')

addArchRGenome('hg38')

AML.input.files = createArrowFiles(inputFiles = c('./data/10026/AML1010_1/fragments.tsv.gz',
                                                  './data/10026/AML1010_23/fragments.tsv.gz',
                                                  './data/10026/AML1010_45/fragments.tsv.gz',
                                                  './data/10026/AML1012_12/fragments.tsv.gz',
                                                  './data/10026/AML1012_34/fragments.tsv.gz',
                                                  './data/10026/AML1026_12/fragments.tsv.gz',
                                                  './data/10026/AML1026_34/fragments.tsv.gz'),
                                     sampleNames = c('AML1010_1','AML1010_23','AML1010_45','AML1012_12', 'AML1012_34','AML1026_12','AML1026_34'),
                                     minTSS = 4, minFrags = 1000, addTileMat = T, addGeneScoreMat = T, threads = 12)

doubScores <- addDoubletScores(input = AML.input.files, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

# create ArchR Project
AML.asap = ArchRProject(ArrowFiles = AML.input.files, outputDirectory = './data/10026/AML.asap', copyArrows = F)

# filter doublets
AML.asap = filterDoublets(AML.asap)

# basic clustering to identify major cell types
AML.asap = addIterativeLSI(ArchRProj = AML.asap, useMatrix = 'TileMatrix', name = 'IterativeLSI', 
                                 iterations = 4, clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
                                 varFeatures = 25000, dimsToUse = 1:30)
AML.asap = addClusters(input = AML.asap, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 1)
set.seed(1987)
AML.asap = addUMAP(ArchRProj = AML.asap, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
AML.asap = addImputeWeights(AML.asap)

saveArchRProject(AML.asap)
AML.asap = loadArchRProject('./data/10026/AML.asap/')

### demultiplex hashing information
source('./R/import_kite_counts.R')

AML.asap$Sample.hash = 'none'
AML.asap$Sample.hash[which(AML.asap$Sample == 'AML1010_1')] = 'AML1010_1'

samples = c('AML1010_23', 'AML1010_45', 'AML1012_12', 'AML1012_34', 'AML1026_12', 'AML1026_34')
samples.1 = c('AML1010_2', 'AML1010_4', 'AML1012_1', 'AML1012_3', 'AML1026_1', 'AML1026_3')
samples.2 = c('AML1010_3', 'AML1010_5', 'AML1012_2', 'AML1012_4', 'AML1026_2', 'AML1026_4')
HTA12.limits = c(1,1,0.5,0.5,0.6,0.8)
HTA34.limits = c(0.5, 0.7, 1, 0.5, 0.6, 0.8)

for (i in seq(1,6)) {
  cellnames = stringr::str_split_fixed(AML.asap$cellNames[which(AML.asap$Sample == samples[i])], pattern = '#', n=2)[,2]
  cellnames = gsub(pattern = '-1', replacement = '', cellnames)
  
  TSA.mat=import_kite_counts(paste0('./data/10026/',samples[i],'.asap/TSA/featurecounts/'))
  TSA.mat = TSA.mat[,which(colnames(TSA.mat) %in% cellnames)]
  
  TSA.so = CreateSeuratObject(TSA.mat, assay = 'HTO')
  TSA.so <- NormalizeData(TSA.so, assay = "HTO", normalization.method = "CLR")
  TSA.so <- ScaleData(TSA.so, assay = "HTO")
  
  boo=as.data.frame(t(GetAssayData(TSA.so)))
  boo$HTA12 = rowMeans(boo[,c('HTA1','HTA2')])
  boo$HTA34 = rowMeans(boo[,c('HTA3','HTA4')])
  p=ggplot(data=boo, aes(x=HTA12, y=HTA34)) + 
    ggrastr::rasterize(ggpointdensity::geom_pointdensity(size=0.5), dpi=600) + 
    scale_color_viridis_c() +
    geom_vline(xintercept = HTA12.limits[i], color='darkred') + 
    geom_hline(yintercept = HTA34.limits[i], color='darkred') +
    theme_classic() + 
    theme(legend.position = 'none',
          axis.title = element_text('Arial', size=10, color='black'),
          axis.text = element_text('Arial', size=10, color='black'))
  ggsave(filename = paste0('./figure_10026/figures/hashing/20210601_',samples[i],'.svg'), plot = p, width = 1.5, height = 1.5)
  
  cells.1 = which(boo$HTA12 > HTA12.limits[i] & boo$HTA34 < HTA34.limits[i])
  cells.2 = which(boo$HTA12 < HTA12.limits[i] & boo$HTA34 > HTA34.limits[i])
  AML.asap$Sample.hash[which(AML.asap$cellNames %in% paste0(samples[i], '#',rownames(boo)[cells.1],'-1'))] = samples.1[i]
  AML.asap$Sample.hash[which(AML.asap$cellNames %in% paste0(samples[i], '#',rownames(boo)[cells.2],'-1'))] = samples.2[i]
}
saveArchRProject(AML.asap)

cells.mat = cbind(AML.asap$cellNames, AML.asap$Sample.hash)
write.table(cells.mat, file = './data/10026/20210611_cells_by_sample.csv', col.names = c('barcode', 'sample'))

### read TSB data
ADT.data = data.frame()
for (s in sort(unique(AML.asap$Sample))) {
  cellnames = stringr::str_split_fixed(AML.asap$cellNames[which(AML.asap$Sample == s)], pattern = '#', n=2)[,2]
  cellnames = gsub(pattern = '-1', replacement = '', cellnames)
  
  TSB.mat=import_kite_counts(paste0('./data/10026/',s,'.asap/TSB/featurecounts/'))
  TSB.mat=TSB.mat[,which(colnames(TSB.mat) %in% cellnames)]
  
  TSB.so = CreateSeuratObject(TSB.mat, assay = 'ADT')
  TSB.so <- NormalizeData(TSB.so, assay = "ADT", normalization.method = 'CLR', margin = 2)
  TSB.so = ScaleData(TSB.so)
  
  boo = as.data.frame(t(GetAssayData(TSB.so)))
  rownames(boo) = paste0(s,'#',rownames(boo))
  colnames(TSB.mat) = paste0(s,'#',colnames(TSB.mat))
  if (s == 'AML1010_1') {
    ADT.data = boo 
    ADT.mat = TSB.mat
  } else {
    ADT.data = rbind(ADT.data, boo)
    ADT.mat = cbind(ADT.mat, TSB.mat)
  }
}

for (protein in colnames(ADT.data)) {
  AML.asap = addCellColData(AML.asap, data = ADT.data[,protein], cells = paste0(rownames(ADT.data),'-1'), name = protein)
}
saveArchRProject(AML.asap)

saveRDS(ADT.data, file = './data/10026/objects/20210601_ADT.rds')

TSB.so = CreateSeuratObject(ADT.mat, assay = 'ADT')
TSB.so <- NormalizeData(TSB.so, assay = "ADT", normalization.method = 'CLR', margin = 2)
TSB.so = ScaleData(TSB.so)
TSB.so = FindVariableFeatures(TSB.so)
TSB.so = RunPCA(TSB.so)
TSB.so = RunUMAP(TSB.so, dims = 1:20)
TSB.so = FindNeighbors(TSB.so)
TSB.so = FindClusters(TSB.so)
saveRDS(TSB.so, file = './data/10026/objects/20210601_TSB.so.rds')

###
p <- plotEmbedding(
  ArchRProj = AML.asap, 
  name = rownames(TSB.so)
)
plotPDF(plotList = p, name = "ASAP.seq", ArchRProj = AML.asap, addDOC = FALSE, width = 5, height = 5)

p = plotEmbedding(AML.asap, name = 'Sample.hash',
                  pal = c('AML1010_1' = RColorBrewer::brewer.pal(n=6, name = 'Reds')[2],
                          'AML1010_2' = RColorBrewer::brewer.pal(n=6, name = 'Reds')[3],
                          'AML1010_3' = RColorBrewer::brewer.pal(n=6, name = 'Reds')[4],
                          'AML1010_4' = RColorBrewer::brewer.pal(n=6, name = 'Reds')[5],
                          'AML1010_5' = RColorBrewer::brewer.pal(n=6, name = 'Reds')[6],
                          'AML1012_1' = RColorBrewer::brewer.pal(n=5, name = 'Greens')[2],
                          'AML1012_2' = RColorBrewer::brewer.pal(n=5, name = 'Greens')[3],
                          'AML1012_3' = RColorBrewer::brewer.pal(n=5, name = 'Greens')[4],
                          'AML1012_4' = RColorBrewer::brewer.pal(n=5, name = 'Greens')[5],
                          'AML1026_1' = RColorBrewer::brewer.pal(n=5, name = 'Blues')[2],
                          'AML1026_2' = RColorBrewer::brewer.pal(n=5, name = 'Blues')[3],
                          'AML1026_3' = RColorBrewer::brewer.pal(n=5, name = 'Blues')[4],
                          'AML1026_4' = RColorBrewer::brewer.pal(n=5, name = 'Blues')[5],
                          'none' = 'grey'))
plotPDF(plotList = list(p), name = "Samples", ArchRProj = AML.asap, addDOC = FALSE, width = 5, height = 5)

p=plotEmbedding(AML.asap, colorBy = 'GeneScoreMatrix', name = 'CD34')
plotPDF(list(p), name = 'CD34', ArchRProj = AML.asap, addDOC = F)

AML.asap <- addGroupCoverages(AML.asap, groupBy = "Clusters", threads = 12)
AML.asap = addReproduciblePeakSet(AML.asap, groupBy = 'Clusters', threads = 12)
AML.asap = addPeakMatrix(AML.asap, threads = 12)
AML.asap <- addMotifAnnotations(AML.asap, motifSet = "cisbp", name = "Motif")
AML.asap <- addBgdPeaks(AML.asap)
AML.asap <- addDeviationsMatrix(AML.asap, peakAnnotation = "Motif", threads = 12, force = T)
saveArchRProject(AML.asap)

# TF motifs
plotVarDev <- getVarDeviations(AML.asap, name = "MotifMatrix", plot = TRUE)
motif.deviations = plotVarDev[["data"]]
motifs = paste0('z:',motif.deviations$name[which(motif.deviations$combinedVars > 5)])
motif.plots = plotEmbedding(AML.asap, colorBy = 'MotifMatrix', name = motifs, embedding = 'UMAP', threads = 1)
plotPDF(motif.plots, name = 'ChromVAR-deviation-scores', ArchRProj = AML.asap, addDOC = F)

p=plotEmbedding(AML.asap, colorBy = 'MotifMatrix', name = 'z:TCF7_750')
plotPDF(list(p), name = 'TCF7_750', ArchRProj = AML.asap, addDOC = F)
p=plotEmbedding(AML.asap, colorBy = 'MotifMatrix', name = 'z:EOMES_788')
plotPDF(list(p), name = 'EOMES_788', ArchRProj = AML.asap, addDOC = F)
p=plotEmbedding(AML.asap, colorBy = 'MotifMatrix', name = 'z:GATA1_383')
plotPDF(list(p), name = 'GATA1_383', ArchRProj = AML.asap, addDOC = F)
p=plotEmbedding(AML.asap, colorBy = 'MotifMatrix', name = 'z:SPI1_322')
plotPDF(list(p), name = 'SPI1_322', ArchRProj = AML.asap, addDOC = F)

# differential TF motifs
markersTF <- getMarkerFeatures(
  ArchRProj = AML.asap, 
  useMatrix = "MotifMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersTF, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 0.5")