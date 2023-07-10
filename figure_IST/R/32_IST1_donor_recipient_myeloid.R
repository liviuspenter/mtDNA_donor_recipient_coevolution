# comparison of donor- and recipient derived myeloid cells in bone marrow of patient IST1

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

donor.cells = rownames(IST.asap.mito)[which(IST.asap.mito$Sample == 'IST2_1' & IST.asap.mito$manual.cluster == 'Myeloid' & IST.asap.mito$individual == 'donor')]
recipient.cells = rownames(IST.asap.mito)[which(IST.asap.mito$Sample == 'IST2_1' & IST.asap.mito$manual.cluster == 'Myeloid' & IST.asap.mito$individual == 'recipient')]

donor.cells.2 = rownames(IST.asap.mito)[which(IST.asap.mito$Sample == 'IST2_2' & IST.asap.mito$manual.cluster == 'Myeloid' & IST.asap.mito$individual == 'donor')]
recipient.cells.2 = rownames(IST.asap.mito)[which(IST.asap.mito$Sample == 'IST2_2' & IST.asap.mito$manual.cluster == 'Myeloid' & IST.asap.mito$individual == 'recipient')]

# downsample donor.cells.2
donor.cells.2 = donor.cells.2[sample(length(donor.cells.2), 100)]

TSB.so = readRDS('./data/IST/objects/20220118_IST_mito_TSB.so.rds')

TSB.so$group = 'none'
TSB.so$group[donor.cells] = 'donor'
TSB.so$group[recipient.cells] = 'recipient'
TSB.so$group[donor.cells.2] = 'donor.2'
TSB.so$group[recipient.cells.2] = 'recipient.2'

VlnPlot(subset(TSB.so, group %in% c('donor', 'recipient', 'donor.2')), 
        group.by = 'group', pt.size = 0.5, slot = 'scale.data',
        features = 'CD33', cols = c('recipient' = 'purple', 'donor' = 'orange', 'donor.2' = 'orange')) +
  scale_x_discrete(limits = c('recipient', 'donor', 'donor.2'),
                   labels = c('recipient pre-IST','donor pre-IST', 'donor post-IST')) + 
  NoLegend() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        plot.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_IST/figures/IST1/plots/20230419_CD33_donor_recipient.svg', width = 1.3, height = 2.5)

VlnPlot(subset(TSB.so, group %in% c('donor', 'recipient', 'donor.2')), 
        group.by = 'group', pt.size = 0.5, slot = 'scale.data',
        features = 'CD16', cols = c('recipient' = 'purple', 'donor' = 'orange', 'donor.2' = 'orange')) +
  scale_x_discrete(limits = c('recipient', 'donor', 'donor.2'),
                   labels = c('recipient pre-IST','donor pre-IST', 'donor post-IST')) + 
  NoLegend() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        plot.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_IST/figures/IST1/plots/20230419_CD16_donor_recipient.svg', width = 1.3, height = 2.5)

VlnPlot(subset(TSB.so, group %in% c('donor', 'recipient', 'donor.2')), 
        group.by = 'group', pt.size = 0.5, slot = 'scale.data',
        features = 'CD14', cols = c('recipient' = 'purple', 'donor' = 'orange', 'donor.2' = 'orange')) +
  scale_x_discrete(limits = c('recipient', 'donor', 'donor.2'),
                   labels = c('recipient pre-IST','donor pre-IST', 'donor post-IST')) + 
  NoLegend() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        plot.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_IST/figures/IST1/plots/20230419_CD14_donor_recipient.svg', width = 1.3, height = 2.5)

VlnPlot(subset(TSB.so, group %in% c('donor', 'recipient', 'donor.2')), 
        group.by = 'group', pt.size = 0.5, slot = 'scale.data',
        features = 'CD11c', cols = c('recipient' = 'purple', 'donor' = 'orange', 'donor.2' = 'orange')) +
  scale_x_discrete(limits = c('recipient', 'donor', 'donor.2'),
                   labels = c('recipient pre-IST','donor pre-IST', 'donor post-IST')) + 
  NoLegend() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        plot.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_IST/figures/IST1/plots/20230419_CD11c_donor_recipient.svg', width = 1.3, height = 2.5)

# differential protein expression
FindMarkers(TSB.so, group.by = 'group', ident.1 = 'donor', ident.2 = 'recipient', logfc.threshold = 0.01, slot = 'scale.data')

protein.mat = GetAssayData(TSB.so, slot = 'scale.data')
protein.mat = protein.mat[, c(recipient.cells, donor.cells, donor.cells.2)]

col_fun = circlize::colorRamp2(breaks = seq(-2,2,4/8), colors = BuenColors::jdb_palette(name = 'brewer_yes',n=9))
ha = columnAnnotation(individual = c(rep('recipient', length(recipient.cells)),
                                     rep('donor', length(donor.cells)),
                                     rep('donor.2', length(donor.cells.2))),
                      col = list('individual' = c('donor' = 'orange', 'donor.2' = 'orange','recipient' = 'purple')),
                      simple_anno_size=unit(5, 'pt'), border=T, annotation_name_gp = gpar(fontsize=0))

TSB.heatmap = Heatmap(protein.mat[c('CD117','CD38','CD33','CD14','CD16', 'CD11c'), 
                                  c(recipient.cells, donor.cells, donor.cells.2)], 
                      show_row_names = T, show_column_names = F, 
                      show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F, 
                      column_split = factor(c(rep('recipient',length(recipient.cells)), 
                                              rep('donor', length(donor.cells)),
                                              rep('donor.2', length(donor.cells.2))),
                                            levels = c('recipient', 'donor', 'donor.2')), top_annotation = ha, 
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
      log2(mean(x[c(donor.cells, donor.cells.2)])) - log2(mean(x[recipient.cells])) 
    }),
    p.value = apply(GS_matrix, MARGIN = 1, FUN = function(x) { 
      t.test(x[c(donor.cells, donor.cells.2)], x[recipient.cells])$p.value 
    })
  )  
volcano.data.GS$p.adj = p.adjust(volcano.data.GS$p.value)
volcano.data.GS$gene = rownames(volcano.data.GS)

volcano.data.GS = volcano.data.GS[-which(volcano.data.GS$log2fc %in% c('-Inf', 'Inf', "NaN")),]

diff.GS = rownames(volcano.data.GS)[which(abs(volcano.data.GS$log2fc) > 1.5 & -log10(volcano.data.GS$p.adj) > 5)]

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
ggsave('./figure_IST/figures/IST1/plots/20230417_IST1_Myeloid_donor_recipient_GS.svg', width = 2, height = 2)

col_fun = circlize::colorRamp2(breaks = seq(0,2,2/8), colors = BuenColors::jdb_palette(name = 'solar_extra',n=9))

GS.heatmap = Heatmap(GS_matrix[c(diff.GS), c(recipient.cells, donor.cells, donor.cells.2)], 
                     show_row_names = T, show_column_names = F, 
                     show_row_dend = F, show_column_dend = F, cluster_rows = T, cluster_columns = F, 
                     column_split = factor(c(rep('recipient',length(recipient.cells)), 
                                             rep('donor', length(donor.cells)),
                                             rep('donor', length(donor.cells.2))),
                                           levels = c('recipient', 'donor', 'donor.2')), 
                     col=col_fun, border=T, row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                     use_raster = T, raster_quality = 10) 

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
      mean(x[c(donor.cells, donor.cells.2)]) - mean(x[recipient.cells]) 
    }),
    p.value = apply(zscore_part, MARGIN = 1, FUN = function(x) { 
      wilcox.test(x[c(donor.cells, donor.cells.2)], x[recipient.cells])$p.value 
    })
  )  
volcano.data.TF$p.adj = p.adjust(volcano.data.TF$p.value)
volcano.data.TF$gene = rownames(volcano.data.TF)

diff.TF = rownames(volcano.data.TF)[which(abs(volcano.data.TF$log2fc) > 1 & -log10(volcano.data.TF$p.adj) > 3)]

ggplot(volcano.data.TF, aes(x=log2fc, y=-log10(p.adj))) + 
  ggrastr::rasterize(geom_point(size=0.5), dpi=600) +
  geom_point(data=volcano.data.TF[diff.TF,], aes(x=log2fc, y=-log10(p.adj)), size=0.5, color='magenta') + 
  geom_label_repel(data=volcano.data.TF[diff.TF,], aes(x=log2fc, y=-log10(p.adj), label=gene), 
                   label.size = 0, size=2, max.overlaps = 15) + 
  scale_x_continuous('Log2FC', limits = c(-2, 2)) + 
  scale_y_continuous('-log10(FDR)') + 
  theme_classic() + 
  theme(axis.title = element_text('Arial', size=8, color='black'),
        axis.text = element_text('Arial', size=8, color='black'))
ggsave('./figure_IST/figures/IST1/plots/20230417_IST1_Myeloid_donor_recipient_TF.svg', width = 2, height = 2)

col_fun = circlize::colorRamp2(breaks = seq(-2,2,4/8), colors = blueYellow)
TF.heatmap = Heatmap(zscore_part[diff.TF, c(recipient.cells, donor.cells, donor.cells.2)], 
                     show_row_names = T, show_column_names = F, 
                     show_row_dend = F, show_column_dend = F, cluster_rows = T, cluster_columns = F, 
                     column_split = factor(c(rep('recipient',length(recipient.cells)), 
                                             rep('donor', length(donor.cells)),
                                             rep('donor.2', length(donor.cells.2))),
                                           levels = c('recipient', 'donor', 'donor.2')), 
                     col=col_fun, border=T, row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                     use_raster = T, raster_quality = 10)


# mtDNA
combined.mutation.frequencies = readRDS('./data/IST/mtDNA/20220110_IST1_2_combined_mutation_frequencies.rds')
germline.variants = read.csv2(file='./data/IST/objects/20220117_IST_germline_variants.csv')

volcano.data.mtDNA = 
  data.frame(
    log2fc = apply(combined.mutation.frequencies, MARGIN = 1, FUN = function(x) { 
      mean(x[c(donor.cells, donor.cells.2)]) - mean(x[recipient.cells]) 
    }),
    p.value = apply(combined.mutation.frequencies, MARGIN = 1, FUN = function(x) { 
      wilcox.test(x[c(donor.cells, donor.cells.2)], x[recipient.cells])$p.value 
    })
  )  
volcano.data.mtDNA$p.adj = p.adjust(volcano.data.mtDNA$p.value)
volcano.data.mtDNA$variant = rownames(volcano.data.mtDNA)

mutations.1 = volcano.data.mtDNA$variant[which(volcano.data.mtDNA$p.adj < 0.05)]
mutations.1 = setdiff(mutations.1, germline.variants$variant[which(germline.variants$sample == 'IST1')])
boo = combined.mutation.frequencies[mutations.1, c(recipient.cells, donor.cells, donor.cells.2)]

mutations.1 = c('3919T>C', '5458T>C', '7457G>A', '10776T>C', '1793G>A', '2623A>G', '8995G>A','13326T>C','13785C>T')
boo = combined.mutation.frequencies[mutations.1, c(recipient.cells, donor.cells, donor.cells.2)]

col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), 
                               colors = BuenColors::jdb_palette(name = 'solar_rojos'))
mtDNA.heatmap = Heatmap(boo, show_row_names = T, show_column_names = F, 
                        show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F, 
                        column_split = factor(c(rep('recipient',length(recipient.cells)), 
                                                rep('donor', length(donor.cells)),
                                                rep('donor.2', length(donor.cells.2))),
                                              levels = c('recipient', 'donor', 'donor.2')), 
                        col=col_fun, border=T, row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                        use_raster = T, raster_quality = 10)


# including some maternal variants and focused on relevant mtDNA mutationss
boo = combined.mutation.frequencies[c('5458T>C', '7457G>A', 
                                      '1793G>A', '2623A>G', '8995G>A', '13326T>C', '13785C>T',
                                      '16093T>C','16224T>C', '16311T>C', 
                                      '16234C>T','16290C>T','16324T>C'), 
                                    c(recipient.cells, donor.cells, donor.cells.2)]

mtDNA.heatmap.2 = Heatmap(boo, show_row_names = T, show_column_names = F, 
                          show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F, 
                          column_split = factor(c(rep('recipient',length(recipient.cells)), 
                                                  rep('donor', length(donor.cells)),
                                                  rep('donor.2', length(donor.cells.2))),
                                                levels = c('recipient', 'donor', 'donor.2')), 
                          col=col_fun, border=T, row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                          use_raster = T, raster_quality = 10)



svglite::svglite('./figure_IST/figures/IST1/heatmaps/20230419_IST1_donor_recipient_Myeloid.svg', width = 4, height = 3.5)
draw(TSB.heatmap %v% GS.heatmap %v% TF.heatmap %v% mtDNA.heatmap)
dev.off()

svglite::svglite('./figure_IST/figures/IST1/heatmaps/20230623_IST1_donor_recipient_Myeloid.svg', width = 4, height = 3.5)
draw(TSB.heatmap %v% GS.heatmap %v% TF.heatmap %v% mtDNA.heatmap.2)
dev.off()

# browser tracks
IST.asap.mito$boo = 'none'
IST.asap.mito$boo[which(rownames(IST.asap.mito) %in% recipient.cells)] = 'recipient'
IST.asap.mito$boo[which(rownames(IST.asap.mito) %in% donor.cells)] = 'donor'
IST.asap.mito$boo[which(rownames(IST.asap.mito) %in% donor.cells.2)] = 'donor.2'

p=plotBrowserTrack(IST.asap.mito, groupBy = 'boo', 
                   useGroups = c('recipient', 'donor', 'donor.2'), geneSymbol = 'IL1B',
                   pal = c('donor' = 'orange', 'donor.2' = 'orange', 'recipient' = 'purple'), 
                   facetbaseSize = 8)
svglite::svglite('./figure_IST/figures/IST1/browser_plots/20230425_IST1_monocytes_IL1B.svg', width = 3, height = 2)
grid.draw(p$IL1B)
dev.off()

# visualize monocytes donor/recipent on UMAP
df = getEmbedding(IST.asap.mito)
colnames(df) = c('UMAP1', 'UMAP2')
p=ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(size=0.5, color='grey90') +
  geom_point(data=df[c(donor.cells, donor.cells.2),], aes(x=UMAP1, y=UMAP2), color='orange', size=2) +
  geom_point(data=df[recipient.cells,], aes(x=UMAP1, y=UMAP2), color='purple', size=2) +
  theme_classic() +
  NoAxes()
ggsave('./figure_IST/figures/IST1/UMAPs/20230425_IST1_monocytes_donor_recipient.png', width = 4, height = 4, dpi = 600, plot = p)



