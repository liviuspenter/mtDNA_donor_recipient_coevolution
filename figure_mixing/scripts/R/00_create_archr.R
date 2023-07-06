setwd('/Users/shaka87/dfci/asap_seq/')

# create ArchR object of CLL4 (CLL_relapse1_1) and CLL5 (CLL_relapse3_1) 
# for extraction of cell barcodes for generation of synthetic bam files

library(ArchR)
library(ggplot2)
library(parallel)
library(Seurat)
library(Signac)

# create ArchR object and put into ArchRProjects storage
setwd('/Users/shaka87/ArchRProjects/CLL_relapse/')

addArchRGenome('hg38')

# path needs to be changed
input.files = createArrowFiles(inputFiles = c('/Users/shaka87/dfci/mito_atac/data/20200730/CLL_relapse1_1/fragments.tsv.gz',
                                              '/Users/shaka87/dfci/mito_atac/data/20200730/CLL_relapse3_1/fragments.tsv.gz'),
                               sampleNames = c('CLL_relapse1_1','CLL_relapse3_1'),
                               minTSS = 4, minFrags = 1000, addTileMat = T, addGeneScoreMat = T)
doubScores <- addDoubletScores(input = input.files, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

CLL.mix = ArchRProject(ArrowFiles = input.files, outputDirectory = './CLL.mix', copyArrows = F)
CLL.mix = filterDoublets(CLL.mix)

CLL.mix = addIterativeLSI(ArchRProj = CLL.mix, useMatrix = 'TileMatrix', name = 'IterativeLSI', 
                                iterations = 4, clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
                                varFeatures = 25000, dimsToUse = 1:30)
CLL.mix = addClusters(input = CLL.mix, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 1)
CLL.mix = addUMAP(ArchRProj = CLL.mix, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 1, metric = 'cosine')
CLL.mix = addImputeWeights(CLL.mix)
CLL.mix$manual.cluster = 'CLL1'
CLL.mix$manual.cluster[which(CLL.mix$Clusters %in% c('C10', 'C11', 'C12'))] = 'CLL3'
CLL.mix$manual.cluster[which(CLL.mix$Clusters %in% c('C1'))] = 'Mono'
CLL.mix$manual.cluster[which(CLL.mix$Clusters %in% c('C2'))] = 'T cell'

saveArchRProject(CLL.mix)
