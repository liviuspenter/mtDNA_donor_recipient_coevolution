# comparison of donor- and recipient derived T cells in peripheral blood of patient IST4 (IST5_1, IST5_2)

library(ArchR)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(gplots)
library(grid)
library(parallel)
library(Signac)
library(Seurat)
library(BuenColors)
library(dplyr)

IST.asap.mito = loadArchRProject('./data/IST/IST.asap.mito/')

donor.cells = rownames(IST.asap.mito)[which(IST.asap.mito$Sample == 'IST5_1' & IST.asap.mito$manual.cluster == 'T' & IST.asap.mito$individual == 'donor')]
recipient.cells = rownames(IST.asap.mito)[which(IST.asap.mito$Sample == 'IST5_1' & IST.asap.mito$manual.cluster == 'T' & IST.asap.mito$individual == 'recipient')]

donor.cells.2 = rownames(IST.asap.mito)[which(IST.asap.mito$Sample == 'IST5_2' & IST.asap.mito$manual.cluster == 'T' & IST.asap.mito$individual == 'donor')]
recipient.cells.2 = rownames(IST.asap.mito)[which(IST.asap.mito$Sample == 'IST5_2' & IST.asap.mito$manual.cluster == 'T' & IST.asap.mito$individual == 'recipient')]

TSB.so = readRDS('./data/IST/objects/20220118_IST_mito_TSB.so.rds')

TSB.so$group = 'none'
TSB.so$group[donor.cells] = 'donor'
TSB.so$group[recipient.cells] = 'recipient'
TSB.so$group[donor.cells.2] = 'donor.2'
TSB.so$group[recipient.cells.2] = 'recipient.2'

# order T cells by CD4
recipient.cells.2 = names(sort(GetAssayData(TSB.so, slot = 'scale.data')['CD4', recipient.cells.2], decreasing = T))
donor.cells.2 = names(sort(GetAssayData(TSB.so, slot = 'scale.data')['CD4', donor.cells.2], decreasing = T))

FeatureScatter(TSB.so, group.by = 'group', pt.size = 0.5, slot = 'scale.data',
               cells = c('recipient' = recipient.cells, 'donor' = donor.cells), 
               feature1 = 'CD4', feature2 = 'CD8', cols = c('recipient' = 'purple', 'donor' = 'orange')) + 
  NoLegend() + 
  theme(plot.title = element_blank(),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_IST/figures/IST4/plots/20230417_IST4_CD4_CD8_donor_recipient.svg', width = 1.5, height = 1.5)

FeatureScatter(TSB.so, group.by = 'group', pt.size = 0.5, slot = 'scale.data',
               cells = c('recipient' = recipient.cells.2, 'donor' = donor.cells.2), 
               feature1 = 'CD4', feature2 = 'CD8', cols = c('recipient.2' = 'purple', 'donor.2' = 'orange')) + 
  NoLegend() + 
  theme(plot.title = element_blank(),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_IST/figures/IST4/plots/20230417_IST4_CD4_CD8_donor_recipient.2.svg', width = 1.5, height = 1.5)

# recipient has more CD4
FindMarkers(TSB.so, group.by = 'group', ident.1 = 'donor.2', ident.2 = 'recipient.2', logfc.threshold = 0.01, slot = 'scale.data')

protein.mat = GetAssayData(TSB.so, slot = 'scale.data')
protein.mat = protein.mat[, c(recipient.cells.2, donor.cells.2)]

col_fun = circlize::colorRamp2(breaks = seq(-2,2,4/8), colors = BuenColors::jdb_palette(name = 'brewer_yes',n=9))
ha = columnAnnotation(individual = c(rep('recipient', length(recipient.cells.2)),
                                     rep('donor', length(donor.cells.2))),
                      col = list('individual' = c('donor' = 'orange', 'recipient' = 'purple')),
                      simple_anno_size=unit(5, 'pt'), border=T, annotation_name_gp = gpar(fontsize=0))

TSB.heatmap = Heatmap(protein.mat[c('CD3','CD4','CD8'), 
                                  c(recipient.cells.2, donor.cells.2)], 
                      show_row_names = T, show_column_names = F, 
                      show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F, 
                      column_split = factor(c(rep('recipient',length(recipient.cells.2)), rep('donor', length(donor.cells.2))),
                                            levels = c('recipient', 'donor')), top_annotation = ha, 
                      column_title_gp = gpar(fontsize=0),
                      col=col_fun, border=T, row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                      use_raster = T, raster_quality = 10)


### GeneScore 
GSMatrix_tmp = getMatrixFromProject(ArchRProj = IST.asap.mito,
                                    useMatrix = "GeneScoreMatrix", useSeqnames = NULL,
                                    verbose = TRUE, binarize = FALSE, threads = 12)

GS_matrix     = as.matrix(assays(GSMatrix_tmp)$GeneScoreMatrix) %>% as.data.frame()
rownames(GS_matrix) = rowData(GSMatrix_tmp)$name

volcano.data.GS = 
  data.frame(
    log2fc = apply(GS_matrix, MARGIN = 1, FUN = function(x) { 
      log2(mean(x[donor.cells.2])) - log2(mean(x[recipient.cells.2])) 
    }),
    p.value = apply(GS_matrix, MARGIN = 1, FUN = function(x) { 
      t.test(x[donor.cells.2], x[recipient.cells.2])$p.value 
    })
  )  
volcano.data.GS$p.adj = p.adjust(volcano.data.GS$p.value)
volcano.data.GS$gene = rownames(volcano.data.GS)

volcano.data.GS = volcano.data.GS[-which(volcano.data.GS$log2fc %in% c('-Inf', 'Inf', "NaN")),]

diff.GS = rownames(volcano.data.GS)[which(abs(volcano.data.GS$log2fc) > 1 & -log10(volcano.data.GS$p.adj) > 2)]

ggplot(volcano.data.GS, aes(x=log2fc, y=-log10(p.adj))) + 
  ggrastr::rasterize(geom_point(size=0.5), dpi=600) +
  geom_point(data=volcano.data.GS[diff.GS,], aes(x=log2fc, y=-log10(p.adj)), size=0.5, color='magenta') + 
  geom_label_repel(data=volcano.data.GS[diff.GS,], aes(x=log2fc, y=-log10(p.adj), label=gene), 
                   label.size = 0, size=2, max.overlaps = 15) + 
  scale_x_continuous('Log2FC') + 
  scale_y_continuous('-log10(FDR)') + 
  theme_classic() + 
  theme(axis.title = element_text('Arial', size=8, color='black'),
        axis.text = element_text('Arial', size=8, color='black'))
ggsave('./figure_IST/figures/IST4/plots/20230417_IST4_Tcells_donor_recipient_GS.svg', width = 2, height = 2)

col_fun = circlize::colorRamp2(breaks = seq(0,2,2/8), colors = BuenColors::jdb_palette(name = 'solar_extra',n=9))

anno = anno_mark(at = which(diff.GS %in% c('PRF1', 'PDCD1', 'GNLY', 'RUNX3')), 
                 labels = diff.GS[which(diff.GS %in% c('PRF1', 'PDCD1', 'GNLY', 'RUNX3'))], which = "row", 
                 labels_gp = gpar(fontsize=8), side = 'left')

GS.heatmap = Heatmap(GS_matrix[c(diff.GS), c(recipient.cells.2, donor.cells.2)], 
                     show_row_names = T, show_column_names = F, 
                     show_row_dend = F, show_column_dend = F, cluster_rows = T, cluster_columns = F, 
                     column_split = factor(c(rep('recipient',length(recipient.cells.2)), rep('donor', length(donor.cells.2))),
                                           levels = c('recipient', 'donor')), 
                     col=col_fun, border=T, row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                     use_raster = T, raster_quality = 10) #+ rowAnnotation(mark = anno)


# MotifMatrix
MotifMatrix_tmp = getMatrixFromProject(ArchRProj = IST.asap.mito,
                                       useMatrix = "MotifMatrix", useSeqnames = NULL,
                                       verbose = TRUE, binarize = FALSE, threads = 12)
zscore_matrix     = as.matrix(assays(MotifMatrix_tmp)$z)          %>% as.data.frame()
deviations_matrix = as.matrix(assays(MotifMatrix_tmp)$deviations) %>% as.data.frame()
VarDeviation = getVarDeviations(IST.asap.mito, name = "MotifMatrix", plot = FALSE) %>% as.data.frame()
Top_TF_target = head(VarDeviation$name, n=150)
zscore_part = zscore_matrix[Top_TF_target,]
zscore_part[zscore_part>2]  = 2
zscore_part[-2>zscore_part] = -2
rownames(zscore_part) = stringr::str_split_fixed(rownames(zscore_part), pattern = '_', n=2)[,1]

volcano.data.TF = 
  data.frame(
    log2fc = apply(zscore_part, MARGIN = 1, FUN = function(x) { 
      mean(x[donor.cells.2]) - mean(x[recipient.cells.2]) 
    }),
    p.value = apply(zscore_part, MARGIN = 1, FUN = function(x) { 
      wilcox.test(x[donor.cells.2], x[recipient.cells.2])$p.value 
    })
  )  
volcano.data.TF$p.adj = p.adjust(volcano.data.TF$p.value)
volcano.data.TF$gene = rownames(volcano.data.TF)

diff.TF = rownames(volcano.data.TF)[which(abs(volcano.data.TF$log2fc) > 1 & -log10(volcano.data.TF$p.adj) > 2)]

ggplot(volcano.data.TF, aes(x=log2fc, y=-log10(p.adj))) + 
  ggrastr::rasterize(geom_point(size=0.5), dpi=600) +
  geom_point(data=volcano.data.TF[diff.TF,], aes(x=log2fc, y=-log10(p.adj)), size=0.5, color='magenta') + 
  geom_label_repel(data=volcano.data.TF[diff.TF,], aes(x=log2fc, y=-log10(p.adj), label=gene), 
                   label.size = 0, size=2, max.overlaps = 15) + 
  scale_x_continuous('Log2FC', limits = c(-1.25, 1.25)) + 
  scale_y_continuous('-log10(FDR)') + 
  theme_classic() + 
  theme(axis.title = element_text('Arial', size=8, color='black'),
        axis.text = element_text('Arial', size=8, color='black'))
ggsave('./figure_IST/figures/IST4/plots/20230417_IST4_Tcells_donor_recipient_TF.svg', width = 2, height = 2)

col_fun = circlize::colorRamp2(breaks = seq(-2,2,4/8), colors = blueYellow)
TF.heatmap = Heatmap(zscore_part[diff.TF, c(recipient.cells.2, donor.cells.2)], 
                     show_row_names = T, show_column_names = F, 
                     show_row_dend = F, show_column_dend = F, cluster_rows = T, cluster_columns = F, 
                     column_split = factor(c(rep('recipient',length(recipient.cells.2)), rep('donor', length(donor.cells.2))),
                                           levels = c('recipient', 'donor')), 
                     col=col_fun, border=T, row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                     use_raster = T, raster_quality = 10)

svglite::svglite('./figure_IST/figures/IST4/heatmaps/20230419_IST4_donor_recipient_Tcells.svg', width = 4, height = 4.5)
draw(TSB.heatmap %v% GS.heatmap %v% TF.heatmap)
dev.off()

VlnPlot(subset(TSB.so, group %in% c('donor.2', 'recipient.2')), pt.size = 0.5, slot = 'scale.data',
               group.by = 'group', features = 'PD1',
               cols = c('donor.2' = 'orange', 'recipient.2' = 'purple')) +
  scale_x_discrete(labels = c('donor', 'recipient')) + 
  NoLegend() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        plot.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_IST/figures/IST4/plots/20230419_PD1_donor_recipient.svg', width = 1.3, height = 2)

boo = data.frame(EOMES = as.numeric(zscore_matrix['EOMES_788', c(donor.cells.2, recipient.cells.2)]),
                 CD4 = GetAssayData(TSB.so, slot = 'scale.data')['CD4', c(donor.cells.2, recipient.cells.2)],
                 CD8 = GetAssayData(TSB.so, slot = 'scale.data')['CD8', c(donor.cells.2, recipient.cells.2)],
                 PD1 = GetAssayData(TSB.so, slot = 'scale.data')['PD1', c(donor.cells.2, recipient.cells.2)])
boo$individual = c(rep('donor', length(donor.cells.2)), rep('recipient', length(recipient.cells.2)))

ggplot(boo, aes(x=EOMES, y=PD1)) + geom_point(aes(color=individual), size=0.5) + 
  scale_color_manual(values = c('donor' = 'orange', 'recipient' = 'purple')) + 
  NoLegend() +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_IST/figures/IST4/plots/20230419_EOMES_PD1_donor_recipient.svg', width = 1.5, height = 1.5)

IST.asap.mito$boo = 'none'
IST.asap.mito$boo[which(rownames(IST.asap.mito) %in% donor.cells.2)] = 'donor'
IST.asap.mito$boo[which(rownames(IST.asap.mito) %in% recipient.cells.2)] = 'recipient'

p=plotBrowserTrack(IST.asap.mito, groupBy = 'boo', useGroups = c('donor', 'recipient'), geneSymbol = 'PRF1')
grid.draw(p$PRF1)

p=plotBrowserTrack(IST.asap.mito, groupBy = 'boo', useGroups = c('donor', 'recipient'), geneSymbol = 'PDCD1', upstream = 10000, downstream = 10000)
grid.draw(p$PDCD1)

p=plotBrowserTrack(IST.asap.mito, groupBy = 'boo', useGroups = c('donor', 'recipient'), geneSymbol = 'RUNX3', 
                   upstream = 100000, downstream = 10000)
grid.draw(p$RUNX3)


df = getEmbedding(IST.asap.mito)
colnames(df) = c('UMAP1', 'UMAP2')
p=ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(size=0.5, color='grey90') +
  geom_point(data=df[c(donor.cells, donor.cells.2),], aes(x=UMAP1, y=UMAP2), color='orange', size=2) +
  geom_point(data=df[c(recipient.cells,recipient.cells.2),], aes(x=UMAP1, y=UMAP2), color='purple', size=2) +
  theme_classic() +
  NoAxes()
ggsave('./figure_IST/figures/IST4/UMAPs/20230425_IST4_Tcells_donor_recipient.png', width = 4, height = 4, dpi = 600, plot = p)

