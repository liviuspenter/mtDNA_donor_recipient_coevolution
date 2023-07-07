setwd('/Users/shaka87/dfci/asap_seq/')

library(ArchR)
library(parallel)
library(Seurat)

MART1.mito = loadArchRProject('/Users/shaka87/ArchRProjects/MART1/MART1.mito/')

ADT.data = readRDS('./data/objects/20211101_ADT.rds')
rownames(ADT.data) = paste0(rownames(ADT.data),'-1')
ADT.data.mito = ADT.data[MART1.mito$cellNames,]
boo = as.data.frame(MART1.mito$donor)
colnames(boo) = 'donor'
rownames(boo) = MART1.mito$cellNames
ADT.data.mito = cbind(ADT.data.mito, boo)

TSB.so = readRDS(file = './data/objects/20211031_MART1_TSB_mito.so.rds')
TSB.so = AddMetaData(TSB.so, metadata = ADT.data.mito)

ADT.data.mito$manual.cluster = TSB.so$manual.cluster[which(colnames(TSB.so) %in% rownames(ADT.data.mito))]

ggplot(ADT.data.mito[which(!ADT.data.mito$manual.cluster %in% c('Mono', 'B cell', 'unspecific')),], 
       aes(x=CD8, y=SAV1, color=donor)) + geom_point() + scale_y_continuous(limits =c(0,10))