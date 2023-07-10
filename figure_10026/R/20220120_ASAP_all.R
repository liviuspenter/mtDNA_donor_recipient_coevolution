library(ArchR)
library(parallel)
library(Rsamtools)
library(gplots)
library(Seurat)

source('/Users/shaka87/dfci/scripts/20200708_filter_doublets_archr_fix.R')
source('/Users/shaka87/dfci/asap_seq/analysis/import_kite_counts.R')
addArchRGenome('hg38')

setwd('/Users/liviuspenter/ArchRProjects/AML.asap/')

ASAP.all.input.files = createArrowFiles(inputFiles = c('/Users/shaka87/dfci/asap_seq/data/AML1010_1/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1010_23/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1010_45/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1012_12/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1012_34/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1026_12/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/AML1026_34/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/IST1_1/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/IST1_2/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/IST2_1/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/IST2_2/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/IST3_1/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/IST3_2/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/IST4_1/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/IST4_2/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/IST5_1/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/IST5_2/fragments.tsv.gz',
                                                  '/Users/shaka87/dfci/asap_seq/data/BM_Mimitou/GSM4732140_Human_BoneMarrow_hg38_fragments.tsv.gz'),
                                   sampleNames = c('AML1010_1','AML1010_23','AML1010_45','AML1012_12', 'AML1012_34','AML1026_12','AML1026_34',
                                                   'IST1_1','IST1_2','IST2_1','IST2_2','IST3_1','IST3_2','IST4_1','IST4_2','IST5_1','IST5_2', 'BM_Mimitou'),
                                   minTSS = 10, minFrags = 1000, addTileMat = T, addGeneScoreMat = T, threads = 1, force=T)

doubScores <- addDoubletScores(input = ASAP.all.input.files, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

# create ArchR Project
AML.asap.all = ArchRProject(ArrowFiles = ASAP.all.input.files, outputDirectory = './AML.ASAP.all', copyArrows = F)

# filter doublets
AML.asap.all = filterDoublets(AML.asap.all)

# basic clustering to identify cells
AML.asap.all = addIterativeLSI(ArchRProj = AML.asap.all, useMatrix = 'TileMatrix', name = 'IterativeLSI', 
                           iterations = 4, clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
                           varFeatures = 25000, dimsToUse = 1:30, threads = 12)
AML.asap.all = addClusters(input = AML.asap.all, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 1)
set.seed(1987)
AML.asap.all = addUMAP(ArchRProj = AML.asap.all, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
AML.asap.all = addImputeWeights(AML.asap.all)
saveArchRProject(AML.asap.all)

setwd('/Users/shaka87/dfci/asap_seq/')

### read TSB data
ADT.data = data.frame()
for (s in sort(unique(AML.asap.all$Sample))) {
  cellnames = stringr::str_split_fixed(AML.asap.all$cellNames[which(AML.asap.all$Sample == s)], pattern = '#', n=2)[,2]
  cellnames = gsub(pattern = '-1', replacement = '', cellnames)
  
  TSB.mat=import_kite_counts2(paste0('./data/',s,'.asap/TSB/featurecounts/'))
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
    common.proteins = intersect(colnames(ADT.data), colnames(boo))
    ADT.data = rbind(ADT.data, boo[,common.proteins])
    ADT.mat = cbind(ADT.mat, TSB.mat[common.proteins,])
  }
}
TSB.so = CreateSeuratObject(ADT.mat, assay = 'ADT')
TSB.so <- NormalizeData(TSB.so, assay = "ADT", normalization.method = 'CLR', margin = 2)
TSB.so = ScaleData(TSB.so)
TSB.so = FindVariableFeatures(TSB.so)
TSB.so = RunPCA(TSB.so)
TSB.so = RunUMAP(TSB.so, dims = 1:10)
TSB.so = FindNeighbors(TSB.so)
TSB.so = FindClusters(TSB.so)
saveRDS(file='./data/objects/20220301_AML_asap_all_TSB.so', TSB.so)

ADT.data[stringr::str_split_fixed(setdiff(AML.asap.all$cellNames, paste0(rownames(ADT.data),'-1')), pattern='-', n=2)[,1],] = 0 
for (protein in colnames(ADT.data)) {
  AML.asap.all = addCellColData(AML.asap.all, data = ADT.data[,protein], cells = paste0(rownames(ADT.data),'-1'), name = protein, force=T)
}

AML.asap.all = subsetArchRProject(AML.asap.all, cells = AML.asap.all$cellNames[which(AML.asap.all$TSSEnrichment > 10)], 
                                  outputDirectory = './AML.ASAP.all.hq/', threads = 12)
AML.asap.all = addClusters(input = AML.asap.all, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 1, force = T)
AML.asap.all = addUMAP(ArchRProj = AML.asap.all, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
AML.asap.all = addImputeWeights(AML.asap.all)
saveArchRProject(AML.asap.all)

AML.asap.all$patient = stringr::str_split_fixed(AML.asap.all$Sample, pattern = '_', n=2)[,1]

AML.asap.all$manual.cluster = 'AML'
AML.asap.all$manual.cluster[which(AML.asap.all$Clusters %in% c('C23', 'C24', 'C25'))] = 'CD8+ Tcells / NK'
AML.asap.all$manual.cluster[which(AML.asap.all$Clusters %in% c('C21', 'C22'))] = 'CD4+ T cells'
AML.asap.all$manual.cluster[which(AML.asap.all$Clusters %in% c('C19'))] = 'B cells'
AML.asap.all$manual.cluster[which(AML.asap.all$Clusters %in% c('C1', 'C11', 'C14', 'C15', 'C16', 'C17'))] = 'Monocytes'
AML.asap.all$manual.cluster[which(AML.asap.all$Clusters %in% c('C18'))] = 'Erythroid'

plotEmbedding(AML.asap.all, name = 'manual.cluster')

AML.asap.all = addGroupCoverages(AML.asap.all)
AML.asap.all = addReproduciblePeakSet(AML.asap.all, groupBy = 'Clusters', threads = 12)
AML.asap.all = addPeakMatrix(AML.asap.all, threads = 12)
AML.asap.all = addMotifAnnotations(AML.asap.all, motifSet = 'cisbp')
AML.asap.all <- addBgdPeaks(AML.asap.all)
