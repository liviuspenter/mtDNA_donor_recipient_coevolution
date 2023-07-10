# comparison of donor- and recipient derived progenitor cells in bone marrow of patient IST1

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

blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")

IST.asap.mito = loadArchRProject('./data/IST/IST.asap.mito/')
combined.mutation.frequencies = readRDS('./data/IST/mtDNA/20220110_IST1_2_combined_mutation_frequencies.rds')
germline.variants = read.csv2(file='./data/IST/objects/20220117_IST_germline_variants.csv')

### comparison of progenitor donor and recipient cells 
# protein
# GeneScore
# TF 
# mtDNA

df = data.frame(bc = IST.asap.mito$cellNames,
                sample = IST.asap.mito$Sample,
                individual = IST.asap.mito$individual,
                manual.cluster = IST.asap.mito$manual.cluster) %>%
  filter(sample %in% c('IST2_1', 'IST2_2'))

donor.cells = c(df$bc[which(df$individual == 'donor' & df$manual.cluster == 'HSC' & df$sample == 'IST2_1')],
                df$bc[which(df$individual == 'donor' & df$manual.cluster == 'HSC' & df$sample == 'IST2_2')])
recipient.cells = df$bc[which(df$individual == 'recipient' & df$manual.cluster == 'HSC')]

# protein
TSB.so = readRDS('./data/IST/objects/20220118_IST_mito_TSB.so.rds')

TSB.so$group = 'none'
TSB.so$group[donor.cells] = 'donor'
TSB.so$group[recipient.cells] = 'recipient'

protein.mat = GetAssayData(TSB.so, slot = 'scale.data')
protein.mat = protein.mat[, c(recipient.cells, donor.cells)]

col_fun = circlize::colorRamp2(breaks = seq(-2,2,4/8), colors = BuenColors::jdb_palette(name = 'brewer_yes',n=9))
ha = columnAnnotation(individual = c(rep('recipient', length(recipient.cells)),
                                     rep('donor', length(donor.cells))),
                      col = list('individual' = c('donor' = 'orange', 'recipient' = 'purple')),
                      simple_anno_size=unit(5, 'pt'), border=T, annotation_name_gp = gpar(fontsize=0))
TSB.heatmap = Heatmap(protein.mat[c('CD117','CD38','CD33','CD14','CD11c', 'CD16'), 
                                  c(recipient.cells, donor.cells)], 
                     show_row_names = T, show_column_names = F, 
                     show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F, 
                     column_split = factor(c(rep('recipient',length(recipient.cells)), rep('donor', length(donor.cells))),
                                           levels = c('recipient', 'donor')), top_annotation = ha, 
                     column_title_gp = gpar(fontsize=0),
                     col=col_fun, border=T, row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                     use_raster = T, raster_quality = 10)

# differential protein expression
FindMarkers(TSB.so, group.by = 'group', ident.1 = 'donor', ident.2 = 'recipient', logfc.threshold = 0.01, slot = 'scale.data')

# GeneScoreMatrix
GSMatrix_tmp = getMatrixFromProject(ArchRProj = IST.asap.mito,
                                    useMatrix = "GeneScoreMatrix", useSeqnames = NULL,
                                    verbose = TRUE, binarize = FALSE, threads = 12)

GS_matrix     = as.matrix(assays(GSMatrix_tmp)$GeneScoreMatrix) %>% as.data.frame()
rownames(GS_matrix) = rowData(GSMatrix_tmp)$name

volcano.data.GS = 
  data.frame(
    log2fc = apply(GS_matrix, MARGIN = 1, FUN = function(x) { 
      log2(mean(x[donor.cells])) - log2(mean(x[recipient.cells])) 
    }),
    p.value = apply(GS_matrix, MARGIN = 1, FUN = function(x) { 
      t.test(x[donor.cells], x[recipient.cells])$p.value 
    })
  )  
volcano.data.GS$p.adj = p.adjust(volcano.data.GS$p.value)
volcano.data.GS$gene = rownames(volcano.data.GS)

volcano.data.GS = volcano.data.GS[-which(volcano.data.GS$log2fc %in% c('-Inf', 'Inf', "NaN")),]

diff.GS = rownames(volcano.data.GS)[which(abs(volcano.data.GS$log2fc) > 1 & -log10(volcano.data.GS$p.adj) > 20)]

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
ggsave('./figure_IST/figures/IST1/plots/20230417_IST1_HSC_donor_recipient_GS.svg', width = 2, height = 2)

col_fun = circlize::colorRamp2(breaks = seq(0,2,2/8), colors = BuenColors::jdb_palette(name = 'solar_extra',n=9))
GS.heatmap = Heatmap(GS_matrix[c(diff.GS), c(recipient.cells, donor.cells)], 
                     show_row_names = T, show_column_names = F, 
                     show_row_dend = F, show_column_dend = F, cluster_rows = T, cluster_columns = F, 
                     column_split = factor(c(rep('recipient',length(recipient.cells)), rep('donor', length(donor.cells))),
                                           levels = c('recipient', 'donor')), 
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
      mean(x[donor.cells]) - mean(x[recipient.cells]) 
    }),
    p.value = apply(zscore_part, MARGIN = 1, FUN = function(x) { 
      t.test(x[donor.cells], x[recipient.cells])$p.value 
    })
  )  
volcano.data.TF$p.adj = p.adjust(volcano.data.TF$p.value)
volcano.data.TF$gene = rownames(volcano.data.TF)

diff.TF = rownames(volcano.data.TF)[which(abs(volcano.data.TF$log2fc) > 0.25 & -log10(volcano.data.TF$p.adj) > 2)]

ggplot(volcano.data.TF, aes(x=log2fc, y=-log10(p.adj))) + 
  ggrastr::rasterize(geom_point(size=0.5), dpi=600) +
  geom_point(data=volcano.data.TF[diff.TF,], aes(x=log2fc, y=-log10(p.adj)), size=0.5, color='magenta') + 
  geom_label_repel(data=volcano.data.TF[diff.TF,], aes(x=log2fc, y=-log10(p.adj), label=gene), 
                   label.size = 0, size=2, max.overlaps = 15) + 
  scale_x_continuous('Log2FC') + 
  scale_y_continuous('-log10(FDR)') + 
  theme_classic() + 
  theme(axis.title = element_text('Arial', size=8, color='black'),
        axis.text = element_text('Arial', size=8, color='black'))
ggsave('./figure_IST/figures/IST1/plots/20230417_IST1_HSC_donor_recipient_TF.svg', width = 2, height = 2)

col_fun = circlize::colorRamp2(breaks = seq(-2,2,4/8), colors = blueYellow)
TF.heatmap = Heatmap(zscore_part[diff.TF, c(recipient.cells, donor.cells)], 
        show_row_names = T, show_column_names = F, 
        show_row_dend = F, show_column_dend = F, cluster_rows = T, cluster_columns = F, 
        column_split = factor(c(rep('recipient',length(recipient.cells)), rep('donor', length(donor.cells))),
                              levels = c('recipient', 'donor')), 
        col=col_fun, border=T, row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
        use_raster = T, raster_quality = 10)

# mtDNA
#mutations.1 = c('3919T>C', '10776T>C', '1793G>A', '2623A>G', '8995G>A','13326T>C','13785C>T')
mutations.1 = c('3919T>C', '5458T>C', '7457G>A', '10776T>C', '1793G>A', '2623A>G', '8995G>A','13326T>C','13785C>T')
boo = combined.mutation.frequencies[mutations.1, c(recipient.cells, donor.cells)]

col_fun = circlize::colorRamp2(breaks = seq(0,0.1,0.1/8), 
                               colors = BuenColors::jdb_palette(name = 'solar_rojos'))
mtDNA.heatmap = Heatmap(boo, show_row_names = T, show_column_names = F, 
                        show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F, 
                        column_split = factor(c(rep('recipient',length(recipient.cells)), rep('donor', length(donor.cells))),
                                              levels = c('recipient', 'donor')), 
                        col=col_fun, border=T, row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                        use_raster = T, raster_quality = 10)

svglite::svglite('./figure_IST/figures/IST1/heatmaps/20230417_IST1_HSC_donor_recipient_heatmap.svg', width = 5, height = 3.5)
draw(TSB.heatmap %v% GS.heatmap %v% TF.heatmap %v% mtDNA.heatmap)
dev.off()