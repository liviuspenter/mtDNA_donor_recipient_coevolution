library(ArchR)
library(Rsamtools)
library(gplots)
library(Seurat)

addArchRGenome('hg38')

AML.asap = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.asap/')

setwd('/Users/liviuspenter/dfci/asap_seq/')

AML.1010.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1010.mito/')
AML.1010.mito$boo = paste0(AML.1010.mito$Sample.hash, '.', AML.1010.mito$manual.clusters)

AML.1012.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1012.mito/')
AML.1012.mito$boo = paste0(AML.1012.mito$Sample.hash, '.', AML.1012.mito$manual.clusters)

AML.1026.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1026.mito/')
AML.1026.mito$boo = paste0(AML.1026.mito$Sample.hash, '.', AML.1026.mito$manual.clusters)

markersGS <- getMarkerFeatures(
  ArchRProj = AML.1012.mito, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "boo",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", useGroups = 'AML1012_2.CD8', bgdGroups = 'AML1012_4.CD8'
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")

p=plotBrowserTrack(AML.1010.mito, groupBy = 'boo', geneSymbol = 'WT1', useGroups = paste0('AML1010_',seq(1,5),'.HSCT'), upstream = 100000, downstream = 100000)
grid.draw(p$WT1)
p=plotBrowserTrack(AML.1012.mito, groupBy = 'boo', geneSymbol = 'WT1', useGroups = paste0('AML1012_',seq(1,4),'.HSCT'), upstream = 100000, downstream = 100000)
grid.draw(p$WT1)
p=plotBrowserTrack(AML.1026.mito, groupBy = 'boo', geneSymbol = 'WT1', useGroups = paste0('AML1026_',seq(1,4),'.HSCT'), upstream = 100000, downstream = 100000)
grid.draw(p$WT1)

p=plotBrowserTrack(AML.1010.mito, groupBy = 'boo', geneSymbol = 'PRAME', useGroups = paste0('AML1010_',seq(1,5),'.HSCT'), upstream = 30000, downstream = 20000)
grid.draw(p$PRAME)
p=plotBrowserTrack(AML.1012.mito, groupBy = 'boo', geneSymbol = 'PRAME', useGroups = paste0('AML1012_',seq(1,4),'.HSCT'), upstream = 30000, downstream = 20000)
grid.draw(p$PRAME)
p=plotBrowserTrack(AML.1026.mito, groupBy = 'boo', geneSymbol = 'PRAME', useGroups = paste0('AML1026_',seq(1,4),'.HSCT'), upstream = 30000, downstream = 20000)
grid.draw(p$PRAME)


