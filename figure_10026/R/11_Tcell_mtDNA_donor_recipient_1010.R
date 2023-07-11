# donor-recipient deconvolution of T cells in AML1010
# used for ASH 2021 abstract (https://doi.org/10.1182/blood-2021-145304)

library(ArchR)
library(gplots)
library(Seurat)
library(ComplexHeatmap)

source('./figure_10026/R/Tcell.colors.R')

AML.Tcell = loadArchRProject('./data/10026/AML.Tcell/')
TSB.Tcell = readRDS('./data/10026/objects/20210603_TSB_Tcell.rds')

### 1010 
s = 'AML1010'
donor.variants = gtools::mixedsort(c('195T>C', '150C>T', '6146A>G', '1811A>G', '10907T>C', '6047A>G', '5999T>C', '14866C>T', '11009T>C', '11467A>G', '9070T>G', '15693T>C', '12308A>G',
                   '12372G>A', '4646T>C', '14620C>T', '15530T>C', '4811A>G', '499G>A', '16356T>C', '11332C>T', '16179C>T'))
recipient.variants = gtools::mixedsort(c('16126T>C', '4216T>C', '10084T>C', '489T>C', '462C>T', '11251A>G', '12612A>G', '15452C>A', '3010G>A', '14798T>C', '13708G>A', 
                       '16069C>T', '16319G>A', '295C>T', '55T>C', '56A>G', '185G>A', '228G>A'))

# read mtDNA data
variants = read.table(paste0('./data/10026/mtDNA/20210429_',s,'_high_confidence_variants.csv'), header = T)
combined.frequencies = readRDS(paste0('./data/10026/mtDNA/20210429_',s,'_combined_mutation_frequencies.rds'))
combined.frequencies = combined.frequencies[-which(rownames(combined.frequencies) == '310T>C'),]

# cells with mtDNA information
cells = AML.Tcell$cellNames[which(grepl(s,AML.Tcell$Sample))]
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

# filter uninformative mtDNA mutations
combined.frequencies = combined.frequencies[-which(rowSums(combined.frequencies) == 0),]
combined.frequencies = combined.frequencies[-which(rowMeans(combined.frequencies) < 0.01),]
#combined.frequencies = combined.frequencies[variants$variant[which(variants$n_cells_conf_detected > 50)],]
#combined.frequencies = combined.frequencies[order(rowVars(as.matrix(combined.frequencies)), decreasing = T)[1:30],]

# order by T cell phenotype
cells.ordered = c()
for (cluster in names(Tcell.colors)) {
  boo = AML.Tcell$cellNames[which(AML.Tcell$manual_clusters == cluster & AML.Tcell$cellNames %in% cells)]
  names(boo) = rep(cluster, length(boo))
  cells.ordered = c(cells.ordered, boo)
}

# score calculated from donor and recipient variants
chimerism.score = as.numeric(colMeans(combined.frequencies[donor.variants,cells.ordered])) - 
  as.numeric(colMeans(combined.frequencies[recipient.variants,cells.ordered]))

chimerism.analysis = ifelse(chimerism.score > 0, 'donor', 'recipient')

# plot heatmap
ha.top = ComplexHeatmap::columnAnnotation(manual_cluster = names(cells.ordered),
                                          chimerism = chimerism.analysis,
                                          show_legend=F,show_annotation_name = T, 
                                          annotation_label = c('manual_cluster' = 'phenotype', 'chimerism' = 'origin'),
                                          annotation_name_gp = gpar(fontsize= 8),
                                          col = list(manual_cluster = Tcell.colors,
                                                     chimerism = c('donor' = 'darkblue', 'recipient' = 'darkred')), border=T)

ha.left = ComplexHeatmap::rowAnnotation(chimerism = c(rep('donor', length(donor.variants)), rep('recipient', length(recipient.variants))), 
                                        show_legend =F, show_annotation_name = F, 
                                        col = list(chimerism = c('donor' = 'darkblue', 'recipient' = 'darkred')), border=T)

col_fun = circlize::colorRamp2(seq(0.01, 1, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))

svglite::svglite('./figure_10026/figures/heatmaps/20210607_AML1010_mtDNA.svg', width = 4, height = 4.5)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[c(donor.variants, recipient.variants),cells.ordered]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = '40 mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = '435 T cells', column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, left_annotation = ha.left, col = col_fun, column_split = chimerism.analysis,
                        row_names_gp = gpar(fontsize = 8), use_raster = T, raster_quality = 10, border=T)
dev.off()

svglite::svglite('./figure_10026/figures/heatmaps/20210607_AML1010_mtDNA.2.svg', width = 3, height = 3)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[c(donor.variants, recipient.variants),cells.ordered]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = F, show_row_names = F, show_heatmap_legend = F, 
                        row_title = '40 mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = '435 T cells', column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, left_annotation = ha.left, col = col_fun, column_split = chimerism.analysis,
                        row_names_gp = gpar(fontsize = 8), use_raster = T, raster_quality = 10, border=T)
dev.off()


# T cells donor versus recipient
donor.cells = cells.ordered[which(chimerism.analysis == 'donor')]
recipient.cells = cells.ordered[which(chimerism.analysis == 'recipient')]
DimPlot(TSB.Tcell, cells.highlight = recipient.cells) + NoLegend()


# chimerism changes
chimerism.df = c(sample = as.character(), cells.sample = as.numeric(), donor.cells = as.numeric(), recipient.cells = as.numeric())
for (s in c('AML1010_1', 'AML1010_2', 'AML1010_3', 'AML1010_4', 'AML1010_5')) {
  cells.sample = AML.Tcell$cellNames[which(AML.Tcell$Sample.hash == s)]
  donor.cells = length(which(chimerism.analysis[which(cells.ordered %in% cells.sample)] == 'donor'))
  recipient.cells = length(which(chimerism.analysis[which(cells.ordered %in% cells.sample)] == 'recipient'))
  chimerism.df = rbind(chimerism.df, data.frame(sample = s, cells.sample = length(cells.sample), 
                                                donor.cells = donor.cells, recipient.cells = recipient.cells))
}
chimerism.df$chimerism = chimerism.df$donor.cells / (chimerism.df$donor.cells + chimerism.df$recipient.cells)
