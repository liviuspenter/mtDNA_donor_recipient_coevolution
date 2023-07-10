library(ArchR)
library(gplots)
library(Seurat)

setwd('/Users/liviuspenter/dfci/asap_seq/')
source('./analysis/Tcell.colors.R')

AML.Tcell = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.Tcell/')

AML.Tcell$boo = 'none'
AML.Tcell$boo[which(AML.Tcell$Sample.hash %in% c('AML1010_1', 'AML1012_1', 'AML1026_1'))] = 'Baseline'
AML.Tcell$boo[which(AML.Tcell$Sample.hash %in% c('AML1010_2', 'AML1012_2', 'AML1026_2'))] = 'EOLN'
AML.Tcell$boo[which(AML.Tcell$Sample.hash %in% c('AML1010_3', 'AML1012_4', 'AML1026_4'))] = 'Ipi'
AML.Tcell$boo2 = paste0(AML.Tcell$boo, '.', AML.Tcell$manual_clusters)

markersGS <- getMarkerFeatures(
  ArchRProj = AML.Tcell, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "boo2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"#, bgdGroups = 'EOLN', useGroups = 'Ipi'
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 0.5")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
  #labelMarkers = c('CD4', 'CD8A', 'PRF1', 'GZMB', 'IL7R', 'TCF7'),
  transpose = F
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")



