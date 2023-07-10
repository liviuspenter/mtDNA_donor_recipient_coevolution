# create ArchR object with all IST samples

library(ArchR)
library(parallel)
library(Rsamtools)
library(gplots)
library(Seurat)

# fix for doublet function - may no longer be necessary in later ArchR versions
source('./R/20200708_filter_doublets_archr_fix.R')

addArchRGenome('hg38')

IST.input.files = createArrowFiles(inputFiles = c('./data/IST/IST1_1/fragments.tsv.gz',
                                                   './data/IST/IST1_2/fragments.tsv.gz',
                                                   './data/IST/IST2_1/fragments.tsv.gz',
                                                   './data/IST/IST2_2/fragments.tsv.gz',
                                                   './data/IST/IST3_1/fragments.tsv.gz',
                                                   './data/IST/IST3_2/fragments.tsv.gz',
                                                   './data/IST/IST4_1/fragments.tsv.gz',
                                                   './data/IST/IST4_2/fragments.tsv.gz',
                                                   './data/IST/IST5_1/fragments.tsv.gz',
                                                   './data/IST/IST5_2/fragments.tsv.gz'),
                                    sampleNames = c('IST1_1','IST1_2','IST2_1','IST2_2','IST3_1','IST3_2','IST4_1','IST4_2','IST5_1','IST5_2'),
                                    minTSS = 4, minFrags = 1000, addTileMat = T, addGeneScoreMat = T, threads = 1)

doubScores <- addDoubletScores(input = IST.input.files, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

# create ArchR Project
IST.asap = ArchRProject(ArrowFiles = IST.input.files, outputDirectory = './data/IST/IST.asap', copyArrows = F)

# filter doublets
IST.asap = filterDoublets(IST.asap)

# basic clustering to identify major cell types
IST.asap = addIterativeLSI(ArchRProj = IST.asap, useMatrix = 'TileMatrix', name = 'IterativeLSI', 
                            iterations = 4, clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
                            varFeatures = 25000, dimsToUse = 1:30)
IST.asap = addClusters(input = IST.asap, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 0.5)
set.seed(1987)
IST.asap = addUMAP(ArchRProj = IST.asap, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
IST.asap = addImputeWeights(IST.asap)

saveArchRProject(IST.asap)
IST.asap = loadArchRProject('./data/IST/IST.asap/')

### read TSB data
source('./R/import_kite_counts.R')
ADT.data = data.frame()
for (s in sort(unique(IST.asap$Sample))) {
  cellnames = stringr::str_split_fixed(IST.asap$cellNames[which(IST.asap$Sample == s)], pattern = '#', n=2)[,2]
  cellnames = gsub(pattern = '-1', replacement = '', cellnames)
  
  TSB.mat=import_kite_counts2(paste0('./data/IST/',s,'.asap/TSB/featurecounts/'))
  TSB.mat=TSB.mat[,which(colnames(TSB.mat) %in% cellnames)]
  
  TSB.so = CreateSeuratObject(TSB.mat, assay = 'ADT')
  TSB.so <- NormalizeData(TSB.so, assay = "ADT", normalization.method = 'CLR', margin = 2)
  TSB.so = ScaleData(TSB.so)
  
  boo = as.data.frame(t(GetAssayData(TSB.so)))
  rownames(boo) = paste0(s,'#',rownames(boo))
  colnames(TSB.mat) = paste0(s,'#',colnames(TSB.mat))
  
  if (s == 'IST1_1') {
    ADT.data = boo 
    ADT.mat = TSB.mat
    TSB.mat.all = TSB.mat
  } else {
    boo[setdiff(colnames(ADT.data), colnames(boo))] = NA
    ADT.data = rbind(ADT.data, boo)
    ADT.mat = merge.sparse(ADT.mat, TSB.mat)
    TSB.mat.all = cbind(TSB.mat.all[intersect(rownames(TSB.mat.all), rownames(TSB.mat)),], TSB.mat[intersect(rownames(TSB.mat.all), rownames(TSB.mat)),])
  }
}
saveRDS(ADT.data, file = './data/IST/objects/20220117_IST_ADT_TSB.rds')

TSB.so = CreateSeuratObject(ADT.mat, assay = 'ADT')
TSB.so <- NormalizeData(TSB.so, assay = "ADT", normalization.method = 'CLR', margin = 2)
TSB.so = ScaleData(TSB.so)
TSB.so = FindVariableFeatures(TSB.so)
TSB.so = RunPCA(TSB.so)
TSB.so = RunUMAP(TSB.so, dims = 1:20)
TSB.so = FindNeighbors(TSB.so)
TSB.so = FindClusters(TSB.so)
saveRDS(TSB.so, file = './data/IST/objects/20220117_IST_TSB.so.rds')

ADT.data[stringr::str_split_fixed(setdiff(IST.asap$cellNames, paste0(rownames(ADT.data),'-1')), pattern='-', n=2)[,1],] = 0 
for (protein in colnames(ADT.data)) {
  IST.asap = addCellColData(IST.asap, data = ADT.data[,protein], cells = paste0(rownames(ADT.data),'-1'), name = protein, force=T)
}
saveArchRProject(IST.asap)
ADT.data.TSB = ADT.data

### read TSA data
ADT.data = data.frame()
for (s in sort(unique(IST.asap$Sample))) {
  cellnames = stringr::str_split_fixed(IST.asap$cellNames[which(IST.asap$Sample == s)], pattern = '#', n=2)[,2]
  cellnames = gsub(pattern = '-1', replacement = '', cellnames)
  
  TSA.mat=import_kite_counts2(paste0('./data/IST/',s,'.asap/TSA/featurecounts/'))
  TSA.mat=TSA.mat[,which(colnames(TSA.mat) %in% cellnames)]
  
  TSA.so = CreateSeuratObject(TSA.mat, assay = 'ADT')
  TSA.so <- NormalizeData(TSA.so, assay = "ADT", normalization.method = 'CLR', margin = 2)
  TSA.so = ScaleData(TSA.so)
  
  boo = as.data.frame(t(GetAssayData(TSA.so)))
  rownames(boo) = paste0(s,'#',rownames(boo))
  colnames(TSA.mat) = paste0(s,'#',colnames(TSA.mat))
  
  if (s == 'IST1_1') {
    ADT.data = boo 
    ADT.mat = TSA.mat
    TSA.mat.all = TSA.mat
  } else {
    boo[setdiff(colnames(ADT.data), colnames(boo))] = 0
    ADT.data[setdiff(colnames(boo), colnames(ADT.data))] = 0
    ADT.data = rbind(ADT.data, boo)
    ADT.mat = merge.sparse(ADT.mat, TSA.mat)
    TSA.mat.all = cbind(TSA.mat.all, TSA.mat[intersect(rownames(TSA.mat.all), rownames(TSA.mat)),])
  }
}
saveRDS(ADT.data, file = './data/IST/objects/20220117_IST_ADT_TSA.rds')

TSA.so = CreateSeuratObject(ADT.mat, assay = 'ADT')
TSA.so <- NormalizeData(TSA.so, assay = "ADT", normalization.method = 'CLR', margin = 2)
TSA.so = ScaleData(TSA.so)
TSA.so = FindVariableFeatures(TSA.so)
TSA.so = RunPCA(TSA.so)
TSA.so = RunUMAP(TSA.so, dims = 1:5)
TSA.so = FindNeighbors(TSA.so, dims = 1:5)
TSA.so = FindClusters(TSA.so)
saveRDS(TSA.so, file = './data/objects/20220117_IST_TSA.so.rds')

ADT.data[stringr::str_split_fixed(setdiff(IST.asap$cellNames, paste0(rownames(ADT.data),'-1')), pattern='-', n=2)[,1],] = 0 
for (protein in colnames(ADT.data)) {
  message(protein)
  IST.asap = addCellColData(IST.asap, data = ADT.data[,protein], cells = paste0(rownames(ADT.data),'-1'), name = protein, force=T)
}
saveArchRProject(IST.asap)