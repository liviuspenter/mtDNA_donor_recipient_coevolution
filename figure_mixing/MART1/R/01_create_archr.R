library(ArchR)
library(ComplexHeatmap)
library(Rsamtools)
library(gplots)
library(Seurat)
library(parallel)

# fix for ArchR - may not be needed in later versions
source('./R/20200708_filter_doublets_archr_fix.R')
# import kite counts from Caleb Lareau
source('./R/import_kite_counts.R')
addArchRGenome('hg38')

MART1.input.files = createArrowFiles(inputFiles = c('./data/mixing/MART1/MART1_1/fragments.tsv.gz',
                                                  './data/mixing/MART1/MART1_2/fragments.tsv.gz',
                                                  './data/mixing/MART1/MART1_3/fragments.tsv.gz'),
                                   sampleNames = c('MART1_1','MART1_2','MART1_3'),
                                   minTSS = 4, minFrags = 1000, addTileMat = T, addGeneScoreMat = T, threads = 1)

doubScores <- addDoubletScores(input = MART1.input.files, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

# create ArchR Project
MART1.asap = ArchRProject(ArrowFiles = MART1.input.files, outputDirectory = './data/mixing/MART1/MART1.asap', copyArrows = F)

# filter doublets
MART1.asap = filterDoublets(MART1.asap)

# basic clustering to identify major cell types
MART1.asap = addIterativeLSI(ArchRProj = MART1.asap, useMatrix = 'TileMatrix', name = 'IterativeLSI', 
                           iterations = 4, clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
                           varFeatures = 25000, dimsToUse = 1:30)
MART1.asap = addClusters(input = MART1.asap, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 1)
set.seed(1987)
MART1.asap = addUMAP(ArchRProj = MART1.asap, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
MART1.asap = addImputeWeights(MART1.asap)

saveArchRProject(MART1.asap)
MART1.asap = loadArchRProject('./data/mixing/MART1/MART1.asap/')

### read TSA data
TSA.mat.all = matrix()
ADT.data = data.frame()
for (s in sort(unique(MART1.asap$Sample))) {
  cellnames = stringr::str_split_fixed(MART1.asap$cellNames[which(MART1.asap$Sample == s)], pattern = '#', n=2)[,2]
  cellnames = gsub(pattern = '-1', replacement = '', cellnames)
  
  TSA.mat=import_kite_counts2(paste0('./data/mixing/MART1/',s,'.asap/TSA/featurecounts/'))
  TSA.mat = TSA.mat[,which(colnames(TSA.mat) %in% cellnames)]
  
  TSA.so = CreateSeuratObject(TSA.mat, assay = 'tetramer')
  TSA.so <- NormalizeData(TSA.so, assay = "tetramer", normalization.method = "CLR")
  TSA.so <- ScaleData(TSA.so, assay = "tetramer")
  
  boo=as.data.frame(t(GetAssayData(TSA.so)))
  rownames(boo) = paste0(s, '#', rownames(boo), '-1')
  colnames(TSA.mat) = paste0(s, '#', colnames(TSA.mat))
  if (s == 'MART1_1') {
    ADT.data = boo 
    TSA.mat.all = TSA.mat
  } else {
    TSA.mat.all = cbind(TSA.mat.all, as.matrix(TSA.mat))
    for (c in setdiff(colnames(ADT.data), colnames(boo))) {
      if (!c %in% colnames(ADT.data)) {
        ADT.data[,c] = NA
      }
      if (!c %in% colnames(boo)) {
        boo[,c] = NA
      }
    }
    ADT.data = rbind(ADT.data, boo)
  }
}

for (cell in MART1.asap$cellNames[which(!MART1.asap$cellNames %in% rownames(ADT.data))]) {
  ADT.data[nrow(ADT.data)+1,] = 0
  rownames(ADT.data)[nrow(ADT.data)] = cell
}
for (protein in colnames(ADT.data)) {
  MART1.asap = addCellColData(MART1.asap, data = ADT.data[,protein], cells = paste0(rownames(ADT.data),''), name = protein)
}
ADT.data.TSA = ADT.data

### read TSB data
ADT.data = data.frame()
for (s in sort(unique(MART1.asap$Sample))) {
  cellnames = stringr::str_split_fixed(MART1.asap$cellNames[which(MART1.asap$Sample == s)], pattern = '#', n=2)[,2]
  cellnames = gsub(pattern = '-1', replacement = '', cellnames)
  
  TSB.mat=import_kite_counts(paste0('./data/mixing/MART1/',s,'.asap/TSB/featurecounts/'))
  TSB.mat=TSB.mat[,which(colnames(TSB.mat) %in% cellnames)]
  
  TSB.so = CreateSeuratObject(TSB.mat, assay = 'ADT')
  TSB.so <- NormalizeData(TSB.so, assay = "ADT", normalization.method = 'CLR', margin = 2)
  TSB.so = ScaleData(TSB.so)
  
  boo = as.data.frame(t(GetAssayData(TSB.so)))
  rownames(boo) = paste0(s,'#',rownames(boo))
  colnames(TSB.mat) = paste0(s,'#',colnames(TSB.mat))
  
  if (s == 'MART1_1') {
    ADT.data = boo 
    ADT.mat = TSB.mat
    TSB.mat.all = TSB.mat
  } else {
    ADT.data = rbind(ADT.data, boo)
    ADT.mat = cbind(ADT.mat, TSB.mat)
    TSB.mat.all = cbind(TSB.mat.all, TSB.mat)
  }
}

for (protein in colnames(ADT.data)) {
  MART1.asap = addCellColData(MART1.asap, data = ADT.data[,protein], cells = paste0(rownames(ADT.data),'-1'), name = protein)
}
saveArchRProject(MART1.asap)

ADT.data = cbind(ADT.data, ADT.data.TSA)
ADT.data$sample = stringr::str_split_fixed(rownames(ADT.data), pattern='#', n=2)[,1]
saveRDS(ADT.data, file = './data/mixing/MART1/objects/20211101_ADT.rds')

TSB.so = CreateSeuratObject(ADT.mat, assay = 'ADT')
TSB.so <- NormalizeData(TSB.so, assay = "ADT", normalization.method = 'CLR', margin = 2)
TSB.so = ScaleData(TSB.so)
TSB.so = FindVariableFeatures(TSB.so)
TSB.so = RunPCA(TSB.so)
TSB.so = RunUMAP(TSB.so, dims = 1:20)
TSB.so = FindNeighbors(TSB.so)
TSB.so = FindClusters(TSB.so)
TSB.so = RenameCells(TSB.so, new.names = paste0(colnames(TSB.so), '-1'))
saveRDS(TSB.so, file = './data/mixing/MART1/objects/20211031_MART1_TSB.so.rds')

### subset cells with available mitochondrial barcoding 
combined.mutation.frequencies = readRDS('./data/mixing/MART1/mtDNA/20211030_MART1_combined_mutation_frequencies.rds')
MART1.mito = subsetArchRProject(ArchRProj = MART1.asap, 
                                   cells = MART1.asap$cellNames[which(MART1.asap$cellNames %in% colnames(combined.mutation.frequencies))],
                                   outputDirectory = './data/mixing/MART1/MART1.mito/', threads = 12, force = T)
MART1.mito = addClusters(input = MART1.mito, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 0.3, force = T)
MART1.mito = addUMAP(ArchRProj = MART1.mito, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
MART1.mito = addImputeWeights(MART1.mito)
saveArchRProject(MART1.mito)
MART1.mito = loadArchRProject('./data/mixing/MART1/MART1.mito/')
  
# plot mitochondrial barcodes
p=plotEmbedding(MART1.mito)
plotPDF(p, name = 'MART1.mito.umap', ArchRProj = MART1.mito, width = 4, height = 4, addDOC = F)
for (mutation in rownames(combined.mutation.frequencies)) {
  MART1.mito$vaf = unlist(combined.mutation.frequencies[mutation, MART1.mito$cellNames])
  MART1.mito$vaf[which(MART1.mito$vaf > 0.1)] = 0.1
  p= plotEmbedding(MART1.mito, name = 'vaf', pal = c('grey','red'), plotAs = 'points', na.rm=T) + 
    ggtitle(mutation) + theme(plot.title = element_text(hjust = 0.5))
  plotPDF(p, name = paste0('MART1.mito.',mutation), ArchRProj = MART1.mito, addDOC = F, width = 4, height = 4)
}

# assign donors - 
# paper annotation: 108 - donor1; 5 - donor2
# donor 108 is enriched in cluster 7 because that is the cluster with CD3 and Tcrb expression
# donor 5 is enriched in cluster 9 because low Tcrb expression

# find relevant maternal mtDNA mutations
boo = MART1.mito$cellNames[which(MART1.mito$Clusters == 'C7')]
mtDNA.108 = names(which(rowMeans(combined.mutation.frequencies[,boo]) > 0.9))
boo = MART1.mito$cellNames[which(MART1.mito$Clusters == 'C9')]
mtDNA.5 = names(which(rowMeans(combined.mutation.frequencies[,boo]) > 0.9))

boo = data.frame(mtDNA.108 = colMeans(combined.mutation.frequencies[mtDNA.108,]),
                 mtDNA.5 = colMeans(combined.mutation.frequencies[mtDNA.5,]))

cells.108 = rownames(boo[which(boo$mtDNA.108 > 0.8 & boo$mtDNA.5 < 0.2),])
cells.5 = rownames(boo[which(boo$mtDNA.108 < 0.2 & boo$mtDNA.5 > 0.8),])
doublets = rownames(boo[which(boo$mtDNA.108 > 0.8 & boo$mtDNA.5 < 0.8),])
unassigned = rownames(boo)[which(!rownames(boo) %in% c(cells.108, cells.5, doublets))]
MART1.mito$donor = NA
MART1.mito$donor[which(MART1.mito$cellNames %in% cells.108)] = '108'
MART1.mito$donor[which(MART1.mito$cellNames %in% cells.5)] = '5'
saveArchRProject(MART1.mito)


p=ggplot() + 
  ggrastr::rasterize(geom_point(data=boo[which(!rownames(boo) %in% c(cells.5, cells.108)),], 
                      aes(x=100*mtDNA.5, y=100*mtDNA.108), size=0.5, color='grey'), dpi=600) + 
  ggrastr::rasterize(geom_point(data=boo[which(rownames(boo) %in% c(cells.5)),], 
             aes(x=100*mtDNA.5, y=100*mtDNA.108), size=0.5, color='purple'), dpi=600) + 
               ggrastr::rasterize(geom_point(data=boo[which(rownames(boo) %in% c(cells.108)),], 
             aes(x=100*mtDNA.5, y=100*mtDNA.108), size=0.5, color='orange'), dpi=600) + 
  scale_x_continuous('% donor 2') +
  scale_y_continuous('% donor 1') +
  geom_vline(xintercept = c(20, 80), color='black') + 
  geom_hline(yintercept = c(20, 80), color='black') + 
  theme_classic() +
  theme(axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/MART1/figures/20230518_mtDNA_donor1_donor2.svg', width = 1.5, height = 1.5, plot = p)


heatmap.cells = c(cells.108[sample(length(cells.108), 500, replace=F)], 
                  cells.5[sample(length(cells.5), 500, replace=F)])
col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = BuenColors::jdb_palette(name = 'solar_rojos', n=9))
ha = columnAnnotation(Donor = c(rep('Donor1', 500), rep('Donor2', 500)), 
                      col = list(Donor = c('Donor1' = 'orange', 'Donor2' = 'blue')), border = T)
svglite::svglite('./figure_mixing/MART1/figures/20211106_MART1_mtDNA_haplotypes.svg', width = 3.5, height = 3.5)
Heatmap(as.matrix(combined.mutation.frequencies[c(gtools::mixedsort(mtDNA.108), gtools::mixedsort(mtDNA.5)), heatmap.cells]), 
        column_split = c(rep('Donor1', 500), rep('Donor2', 500)), cluster_rows = F, show_column_names = F, cluster_columns = F,
        col = col_fun, border=T, row_names_gp = gpar(fontsize=8), column_title_gp = gpar(fontsize=10),
        use_raster = T, raster_quality = 5, top_annotation = ha)
dev.off()

boo = data.frame(sample = c('MART1_1', 'MART1_2', 'MART1_3'),
                 donor = c(length(which(MART1.mito$Sample == 'MART1_1' & MART1.mito$donor == '108')),
                           length(which(MART1.mito$Sample == 'MART1_2' & MART1.mito$donor == '108')),
                           length(which(MART1.mito$Sample == 'MART1_3' & MART1.mito$donor == '108'))),
                 cells = c(length(which(MART1.mito$Sample == 'MART1_1')),
                           length(which(MART1.mito$Sample == 'MART1_2')),
                           length(which(MART1.mito$Sample == 'MART1_3'))))

ggplot(data=boo, aes(x=sample, y=100*donor/cells)) + 
  geom_line(aes(group=1)) + 
  geom_point(aes(size=cells), color='orange') + 
  scale_x_discrete('dilution',labels = c('1:3', '1:30', '1:300')) + 
  scale_y_log10('% donor 108') +
  scale_size_continuous('# all cells',limits= c(4000,6000), range = c(1,3), breaks = c(4000,5000,6000)) + 
  theme_classic() +
  theme(legend.position = 'right',
        axis.title = element_text('Arial', size=10, color='black'),
        #axis.title.x = element_blank(),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/MART1/figures/20211031_MART1_titration.svg', width = 3.5, height = 1.5)  

###
cells = intersect(colnames(TSA.mat.all), colnames(TSB.mat.all))
TSAB.mat = rbind(TSA.mat.all[,cells], TSB.mat.all[,cells])

TSAB.so = CreateSeuratObject(TSAB.mat, assay = 'ADT')
TSAB.so <- NormalizeData(TSAB.so, assay = "ADT", normalization.method = 'CLR', margin = 2)
TSAB.so = ScaleData(TSAB.so)
TSAB.so = FindVariableFeatures(TSAB.so)
TSAB.so = FindNeighbors(TSAB.so)
TSAB.so = FindClusters(TSAB.so, resolution = 2)
TSAB.so = RunPCA(TSAB.so)
TSAB.so = RunUMAP(TSAB.so, dims = 1:10)
saveRDS(file='./data/mixing/MART1/objects/20220225_MART1_TSAB.so', object = TSAB.so)

