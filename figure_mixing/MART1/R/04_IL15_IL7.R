### effect of IL-7/IL-15

library(ArchR)
library(ggplot2)
library(ggrepel)
library(parallel)
library(Seurat)

MART1.mito = loadArchRProject('./data/mixing/MART1/MART1.mito/')

markersGS <- getMarkerFeatures(
  ArchRProj = MART1.asap, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = c('C7', 'C10'),
  bgdGroups = c('C11', 'C9', 'C8', 'C6', 'C5'),
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

MART1.mito = addPeakMatrix(MART1.mito, threads = 12)
MART1.mito = addGroupCoverages(MART1.mito, groupBy = 'Clusters', threads = 1)
MART1.mito = addReproduciblePeakSet(MART1.mito, groupBy = 'Clusters', threads = 12)
MART1.mito = addMotifAnnotations(MART1.mito, motifSet = 'cisbp')
MART1.mito <- addBgdPeaks(MART1.mito)
MART1.mito <- addDeviationsMatrix(MART1.mito, peakAnnotation = "Motif", threads = 1, force=T)
saveArchRProject(MART1.mito)

# TF motifs
plotVarDev <- getVarDeviations(MART1.mito, name = "MotifMatrix", plot = TRUE)
motif.deviations = plotVarDev[["data"]]
motifs = paste0('z:',motif.deviations$name[which(motif.deviations$combinedVars > 4)])
motif.plots = plotEmbedding(MART1.mito, colorBy = 'MotifMatrix', name = motifs, embedding = 'UMAP', threads = 1)
plotPDF(motif.plots, name = 'ChromVAR-deviation-scores', ArchRProj = MART1.mito, addDOC = F)

# identify CD4/8 T cells unstimulated and CD4/8 T cells stimulated
CD4.unstim = MART1.mito$cellNames[which(MART1.mito$Clusters %in% c('C3', 'C6', 'C8', 'C9') & MART1.mito$CD4 > 2.5 & MART1.mito$CD8 < 2)]
CD4.stim = MART1.mito$cellNames[which(MART1.mito$Clusters %in% c('C7') & MART1.mito$CD4 > 2.5 & MART1.mito$CD8 < 2)]
CD8.unstim = MART1.mito$cellNames[which(MART1.mito$Clusters %in% c('C3', 'C6', 'C8', 'C9') & MART1.mito$CD4 < 2 & MART1.mito$CD8 > 2.5)]
CD8.stim = MART1.mito$cellNames[which(MART1.mito$Clusters %in% c('C7') & MART1.mito$CD4 < 2 & MART1.mito$CD8 > 2.5)]
MART1.mito$Tcell.status = 'none'
MART1.mito$Tcell.status[which(MART1.mito$cellNames %in% CD4.unstim)] = 'CD4.unstim'
MART1.mito$Tcell.status[which(MART1.mito$cellNames %in% CD4.stim)] = 'CD4.stim'
MART1.mito$Tcell.status[which(MART1.mito$cellNames %in% CD8.unstim)] = 'CD8.unstim'
MART1.mito$Tcell.status[which(MART1.mito$cellNames %in% CD8.stim)] = 'CD8.stim'

# Differential chromatin analysis
MotifMatrix_tmp = getMatrixFromProject(ArchRProj = MART1.mito,
                                       useMatrix = "MotifMatrix",
                                       useSeqnames = NULL,
                                       verbose = TRUE,
                                       binarize = FALSE,
                                       threads = 12)

zscore_matrix     = as.matrix(assays(MotifMatrix_tmp)$z)          %>% as.data.frame()
deviations_matrix = as.matrix(assays(MotifMatrix_tmp)$deviations) %>% as.data.frame()

VarDeviation = getVarDeviations(MART1.mito, name = "MotifMatrix", plot = FALSE) %>% as.data.frame()
Top_TF_target = head(VarDeviation$name, n=100)
zscore_part = zscore_matrix[Top_TF_target,]
zscore_part[zscore_part>2]  = 2
zscore_part[-2>zscore_part] = -2

# calculate data for volcano plot
clone1 = CD4.stim
clone2 = CD4.unstim
volcano.data = data.frame()
for (TF in rownames(zscore_matrix)) {
  if (length(which(is.na(zscore_matrix[TF,]))) == 0) {
    p.value = wilcox.test(as.numeric(zscore_matrix[TF,clone1]), as.numeric(zscore_matrix[TF,clone2]))$p.value
    diff.accessibility = median(as.numeric(zscore_matrix[TF,clone1])) - median(as.numeric(zscore_matrix[TF,clone2]))
    volcano.data = rbind(volcano.data, data.frame(TF = TF, p.value = p.value, diff.accessibility = diff.accessibility))
  }
}
volcano.data$p.adj = p.adjust(volcano.data$p.value)
volcano.data$p.adj[which(volcano.data$p.adj == 0)] = min(volcano.data$p.adj[which(volcano.data$p.adj != 0)])
write.table(volcano.data, file = './data/mixing/MART1/DGEA/20220228_CD4_stim_unstim.csv')

# calculate data for volcano plot
clone1 = CD8.stim
clone2 = CD8.unstim
volcano.data = data.frame()
for (TF in rownames(zscore_matrix)) {
  if (length(which(is.na(zscore_matrix[TF,]))) == 0) {
    p.value = wilcox.test(as.numeric(zscore_matrix[TF,clone1]), as.numeric(zscore_matrix[TF,clone2]))$p.value
    diff.accessibility = median(as.numeric(zscore_matrix[TF,clone1])) - median(as.numeric(zscore_matrix[TF,clone2]))
    volcano.data = rbind(volcano.data, data.frame(TF = TF, p.value = p.value, diff.accessibility = diff.accessibility))
  }
}
volcano.data$p.adj = p.adjust(volcano.data$p.value)
volcano.data$p.adj[which(volcano.data$p.adj == 0)] = min(volcano.data$p.adj[which(volcano.data$p.adj != 0)])
write.table(volcano.data, file = './data/mixing/MART1/DGEA/20220228_CD8_stim_unstim.csv', sep = ' ', quote = F)

volcano.data.CD4 = read.csv2('./data/mixing/MART1/DGEA/20220228_CD4_stim_unstim.csv', sep = ' ')
volcano.data.CD8 = read.csv2('./data/mixing/MART1/DGEA/20220228_CD8_stim_unstim.csv', sep = ' ')

volcano.data.CD8$diff.accessibility = as.numeric(volcano.data.CD8$diff.accessibility)
volcano.data.CD4$diff.accessibility = as.numeric(volcano.data.CD4$diff.accessibility)
volcano.data.CD4$p.adj = as.numeric(volcano.data.CD4$p.adj)
volcano.data.CD8$p.adj = as.numeric(volcano.data.CD8$p.adj)


volcano.data.CD4$significant = ifelse(abs(volcano.data.CD4$diff.accessibility) > 1.5 & -log10(volcano.data.CD4$p.adj) >50, T, F)
volcano.data.CD8$significant = ifelse(abs(volcano.data.CD8$diff.accessibility) > 1.5 & -log10(volcano.data.CD8$p.adj) >50, T, F)

ggplot() + 
  geom_vline(xintercept = c(-1.5,1.5), color='grey', linetype='dashed') +
  geom_hline(yintercept = c(50), color='grey', linetype='dashed') +
  geom_point(data=volcano.data.CD4[!volcano.data.CD4$significant,], aes(x=diff.accessibility, y=-log10(p.adj)), color='grey') +
  geom_point(data=volcano.data.CD4[volcano.data.CD4$significant,], aes(x=diff.accessibility, y=-log10(p.adj)), color='black') +
  geom_label_repel(data=volcano.data.CD4[volcano.data.CD4$significant,], 
                   aes(x=diff.accessibility, y=-log10(p.adj),label=stringr::str_split_fixed(TF, pattern = '_', n = 2)[,1]), size = 2.5, 
                   label.size = 0, label.padding = 0.1, segment.size = 0.25, segment.colour = 'grey', min.segment.length = 0.1, max.overlaps = 15) + 
  scale_x_continuous('Difference in chromatin accessibility') +
  scale_y_continuous('-log10(FDR)',limits = c(0,150)) + 
  theme_bw() + 
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black', size = 0.25),
        axis.ticks = element_line(colour = 'black', size = 0.25),
        axis.title = element_text(family = 'Arial', size = 10, colour = 'black'),
        axis.text = element_text(family = 'Arial', size = 10, colour = 'black'))
ggsave('./figure_mixing/MART1/figures/20220228_volcano_CD4_stim_unstim.svg', width = 2.5, height = 2)

ggplot() + 
  geom_vline(xintercept = c(-1.5,1.5), color='grey', linetype='dashed') +
  geom_hline(yintercept = c(50), color='grey', linetype='dashed') +
  geom_point(data=volcano.data.CD8[!volcano.data.CD8$significant,], aes(x=diff.accessibility, y=-log10(p.adj)), color='grey') +
  geom_point(data=volcano.data.CD8[volcano.data.CD8$significant,], aes(x=diff.accessibility, y=-log10(p.adj)), color='black') +
  geom_label_repel(data=volcano.data.CD8[volcano.data.CD8$significant,], 
                   aes(x=diff.accessibility, y=-log10(p.adj),label=stringr::str_split_fixed(TF, pattern = '_', n = 2)[,1]), size = 2.5, 
                   label.size = 0, label.padding = 0.1, segment.size = 0.25, segment.colour = 'grey', min.segment.length = 0.1, max.overlaps = 15) + 
  scale_x_continuous('Difference in chromatin accessibility') +
  scale_y_continuous('-log10(FDR)',limits = c(0,300)) + 
  theme_bw() + 
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black', size = 0.25),
        axis.ticks = element_line(colour = 'black', size = 0.25),
        axis.title = element_text(family = 'Arial', size = 10, colour = 'black'),
        axis.text = element_text(family = 'Arial', size = 10, colour = 'black'))
ggsave('./figure_mixing/MART1/figures/20220228_volcano_CD8_stim_unstim.svg', width = 2.5, height = 2)