library(ArchR)
library(gplots)
library(Seurat)
library(ComplexHeatmap)

setwd('/Users/liviuspenter/dfci/asap_seq/')
source('./analysis/Tcell.colors.R')

AML.Tcell = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.Tcell/')

### 1010 
s = 'AML1010'

# read mtDNA data
variants = read.table(paste0('./data/mtDNA/20210429_',s,'_high_confidence_variants.csv'), header = T)
combined.frequencies = readRDS(paste0('./data/mtDNA/20210429_',s,'_combined_mutation_frequencies.rds'))

# cells with mtDNA information
cells = AML.Tcell$cellNames[which(grepl(s,AML.Tcell$Sample))]
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

# filter uninformative mtDNA mutations
combined.frequencies = combined.frequencies[-which(rowSums(combined.frequencies) == 0),]
combined.frequencies = combined.frequencies[-which(rowMeans(combined.frequencies) < 0.001),]
#combined.frequencies = combined.frequencies[variants$variant[which(variants$n_cells_conf_detected > 50)],]
#combined.frequencies = combined.frequencies[order(rowVars(as.matrix(combined.frequencies)), decreasing = T)[1:30],]

# order by T cell phenotype
cells.ordered = c()
for (cluster in names(Tcell.colors)) {
  boo = AML.Tcell$cellNames[which(AML.Tcell$manual_clusters == cluster & AML.Tcell$cellNames %in% cells)]
  names(boo) = rep(cluster, length(boo))
  cells.ordered = c(cells.ordered, boo)
}

# plot heatmap
ha.top = ComplexHeatmap::columnAnnotation(manual_cluster = names(cells.ordered),
                                          show_legend=F,show_annotation_name = F, 
                                          col = list(manual_cluster = Tcell.colors))

col_fun = circlize::colorRamp2(seq(0.01, 1, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))

#svglite::svglite('./figures/heatmaps/BM/20210430_samples.svg', width = 2.5, height = 2.5)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[,cells.ordered]), show_row_dend = F, show_column_dend = T, cluster_columns = T, cluster_rows = T,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = 'Response', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = 'Timepoint', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, col = col_fun, row_names_gp = gpar(fontsize = 10))
#dev.off()


for (s in c('AML1010', 'AML1012', 'AML1026')) {
  #s = 'AML1010'
  
  # read mtDNA data
  variants = read.table(paste0('./data/mtDNA/20210429_',s,'_high_confidence_variants.csv'), header = T)
  combined.frequencies = readRDS(paste0('./data/mtDNA/20210429_',s,'_combined_mutation_frequencies.rds'))
  
  # cells with mtDNA information
  cells = AML.Tcell$cellNames[which(grepl(s,AML.Tcell$Sample))]
  cells = cells[which(cells %in% colnames(combined.frequencies))]
  combined.frequencies = combined.frequencies[,cells]
  
  # filter uninformative mtDNA mutations
  combined.frequencies = combined.frequencies[-which(rowSums(combined.frequencies) == 0),]
  combined.frequencies = combined.frequencies[-which(rowMeans(combined.frequencies) < 0.001),]
  #combined.frequencies = combined.frequencies[variants$variant[which(variants$n_cells_conf_detected > 50)],]
  #combined.frequencies = combined.frequencies[order(rowVars(as.matrix(combined.frequencies)), decreasing = T)[1:30],]
  
  # order by T cell phenotype
  cells.ordered = c()
  for (cluster in names(Tcell.colors)) {
    boo = AML.Tcell$cellNames[which(AML.Tcell$manual_clusters == cluster & AML.Tcell$cellNames %in% cells)]
    names(boo) = rep(cluster, length(boo))
    cells.ordered = c(cells.ordered, boo)
  }
  
  # plot heatmap
  ha.top = ComplexHeatmap::columnAnnotation(manual_cluster = names(cells.ordered),
                                            show_legend=F,show_annotation_name = F, 
                                            col = list(manual_cluster = Tcell.colors))
  
  col_fun = circlize::colorRamp2(seq(0.001, 0.1, length=100), colorRampPalette(c('white', 'red'))(100))
  
  #svglite::svglite('./figures/heatmaps/BM/20210430_samples.svg', width = 2.5, height = 2.5)
  ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[,cells.ordered]), show_row_dend = F, show_column_dend = F, cluster_columns = T, cluster_rows = T,
                          show_column_names = F, show_row_names = F, show_heatmap_legend = F, 
                          row_title = 'Response', row_title_gp = gpar(fontsize=10, color='black'),
                          column_title = 'Timepoint', column_title_gp = gpar(fontsize=10, color='black'),
                          top_annotation = ha.top, col = col_fun, row_names_gp = gpar(fontsize = 10))
  #dev.off()
}