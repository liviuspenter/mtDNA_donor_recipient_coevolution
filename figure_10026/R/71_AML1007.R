library(ArchR)
library(Rsamtools)
library(gplots)

source('./R/20200708_filter_doublets_archr_fix.R')

addArchRGenome('hg38')

AML.input.files = createArrowFiles(inputFiles = c('./data/10026/Pool148_1/fragments.tsv.gz',
                                                  './data/10026/Pool148_2/fragments.tsv.gz',
                                                  './data/10026/Pool148_12/fragments.tsv.gz'),
                                   sampleNames = c('AML1007_1','AML1007_5','AML1007_6'),
                                   minTSS = 4, minFrags = 1000, addTileMat = T, addGeneScoreMat = T, threads = 1)

doubScores <- addDoubletScores(input = AML.input.files, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

# create ArchR Project
AML.1007 = ArchRProject(ArrowFiles = AML.input.files, outputDirectory = './AML.1007', copyArrows = F)

# filter doublets
AML.1007 = filterDoublets(AML.1007)

# basic clustering to identify CLL cells
AML.1007 = addIterativeLSI(ArchRProj = AML.1007, useMatrix = 'TileMatrix', name = 'IterativeLSI', 
                           iterations = 4, clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
                           varFeatures = 25000, dimsToUse = 1:30)
AML.1007 = addClusters(input = AML.1007, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 1)
set.seed(1987)
AML.1007 = addUMAP(ArchRProj = AML.1007, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
AML.1007 = addImputeWeights(AML.1007)

saveArchRProject(AML.1007)

### subset cells with available mitochondrial barcoding 
combined.mutation.frequencies = readRDS('/Users/shaka87/dfci/asap_seq/data/mtDNA/20240516_AML1007_combined_mutation_frequencies.rds')
AML.1007.mito = subsetArchRProject(ArchRProj = AML.1007, 
                                     cells = AML.1007$cellNames[which(AML.1007$cellNames %in% colnames(combined.mutation.frequencies))],
                                     outputDirectory = './AML.1007.mito/', threads = 12)
AML.1007.mito = addClusters(input = AML.1007.mito, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 0.3, force = T)
AML.1007.mito = addUMAP(ArchRProj = AML.1007.mito, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
AML.1007.mito = addImputeWeights(AML.1007.mito)
saveArchRProject(AML.1007.mito)

# plot mitochondrial barcodes
p=plotEmbedding(AML.1007.mito)
plotPDF(p, name = 'AML.1007.mito.umap', ArchRProj = AML.1007.mito, width = 4, height = 4, addDOC = F)
for (mutation in rownames(combined.mutation.frequencies)) {
  AML.1007.mito$vaf = unlist(combined.mutation.frequencies[mutation, AML.1007.mito$cellNames])
  AML.1007.mito$vaf[which(AML.1007.mito$vaf > 0.1)] = 0.1
  p= plotEmbedding(AML.1007.mito, name = 'vaf', pal = c('grey','red'), plotAs = 'points', na.rm=T) + 
    ggtitle(mutation) + theme(plot.title = element_text(hjust = 0.5))
  plotPDF(p, name = paste0('AML.1007.mito.',mutation), ArchRProj = AML.1007.mito, addDOC = F, width = 4, height = 4)
}