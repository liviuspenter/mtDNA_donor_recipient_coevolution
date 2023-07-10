library(ArchR)
library(Rsamtools)
library(gplots)
library(Seurat)

source('/Users/shaka87/dfci/scripts/20200708_filter_doublets_archr_fix.R')

addArchRGenome('hg38')

setwd('/Users/liviuspenter/ArchRProjects/AML.asap/')

AML.input.files = createArrowFiles(inputFiles = c('/Users/shaka87/dfci/asap_seq/data/AML1010_1/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1010_23/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1010_45/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1012_12/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1012_34/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1026_12/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1026_34/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/cimac_9204/10x/BM_35/fragments.tsv.gz'),
                                   sampleNames = c('AML1010_1','AML1010_23','AML1010_45','AML1012_12', 'AML1012_34','AML1026_12','AML1026_34', 'BM35'),
                                   minTSS = 4, minFrags = 1000, addTileMat = T, addGeneScoreMat = T, threads = 12)


# create ArchR Project
AML.asap = ArchRProject(ArrowFiles = AML.input.files, outputDirectory = './AML.asap.with.BM35', copyArrows = F)

# filter doublets
AML.asap = filterDoublets(AML.asap)

# basic clustering to identify  cells
AML.asap = addIterativeLSI(ArchRProj = AML.asap, useMatrix = 'TileMatrix', name = 'IterativeLSI', 
                           iterations = 4, clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
                           varFeatures = 25000, dimsToUse = 1:30)
AML.asap = addClusters(input = AML.asap, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 1)
set.seed(1987)
AML.asap = addUMAP(ArchRProj = AML.asap, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
AML.asap = addImputeWeights(AML.asap)

saveArchRProject(AML.asap)
