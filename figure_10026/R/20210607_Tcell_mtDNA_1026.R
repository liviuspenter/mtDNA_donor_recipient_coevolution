library(ArchR)
library(gplots)
library(Seurat)
library(ComplexHeatmap)

setwd('/Users/liviuspenter/dfci/asap_seq/')
source('./analysis/Tcell.colors.R')

AML.Tcell = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.Tcell/')

### 1026 
s = 'AML1026'
donor.variants = gtools::mixedsort(c('10463T>C', '16519T>C', '15928G>A', '16126T>C', '16153G>A', '73A>G', '15607A>G', '15452C>A', '7028C>T', '4917A>G',
                                     '13368G>A', '4216T>C', '8697G>A', '11812A>G', '14766C>T', '11719G>A', '709G>A', '8269G>A', '14905G>A', '2706A>G',
                                     '11251A>G', '14233A>G', '1888G>A', '9947G>A', '150C>T', '16294C>T', '16296C>T'))
recipient.variants = gtools::mixedsort(c('16304T>C', '456C>T', '8433T>C', '15833C>T', '4336T>C', '9722T>C', '4011C>T', '93A>G'))

# read mtDNA data
variants = read.table(paste0('./data/mtDNA/20210429_',s,'_high_confidence_variants.csv'), header = T)
combined.frequencies = readRDS(paste0('./data/mtDNA/20210429_',s,'_combined_mutation_frequencies.rds'))

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

boo=as.data.frame(chimerism.score)
ggplot(data=boo, aes(x=chimerism.score)) + geom_histogram()

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

svglite::svglite('./figures/heatmaps/20210607_AML1026_mtDNA.svg', width = 4, height = 4.5)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[c(donor.variants, recipient.variants),cells.ordered]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = '35 mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = '796 T cells', column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, left_annotation = ha.left, col = col_fun, column_split = chimerism.analysis,
                        row_names_gp = gpar(fontsize = 8), use_raster = T, raster_quality = 10, border = T)
dev.off()

svglite::svglite('./figures/heatmaps/20210607_AML1026_mtDNA.2.svg', width = 3, height = 3)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[c(donor.variants, recipient.variants),cells.ordered]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = F, show_row_names = F, show_heatmap_legend = F, 
                        row_title = '35 mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = '796 T cells', column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, left_annotation = ha.left, col = col_fun, column_split = chimerism.analysis,
                        row_names_gp = gpar(fontsize = 8), use_raster = T, raster_quality = 10, border = T)
dev.off()


# clusters for donor and recipient
chimerism.df = data.frame(cluster = as.character(), donor=as.numeric(), recipient=as.numeric())
for (cluster in names(Tcell.colors)) {
  donor.cells = length(which(chimerism.analysis[which(names(cells.ordered) == cluster)] == 'donor'))
  recipient.cells = length(which(chimerism.analysis[which(names(cells.ordered) == cluster)] == 'recipient'))
  chimerism.df = rbind(chimerism.df, data.frame(cluster = cluster, donor = donor.cells, recipient = recipient.cells))
}
chimerism.df$donor.freq = chimerism.df$donor / sum(chimerism.df$donor)
chimerism.df$recipient.freq = chimerism.df$recipient / sum(chimerism.df$recipient)
chimerism.df$donor.freq[which(chimerism.df$donor.freq == 0)] = 0.001
chimerism.df$recipient.freq[which(chimerism.df$recipient.freq == 0)] = 0.001

# detect statistically significant differences
chimerism.df$p.value = apply(chimerism.df, 1, function(x) fisher.test(matrix(c(as.numeric(x['recipient']), sum(chimerism.df$recipient), 
                                                        as.numeric(x['donor']), sum(chimerism.df$donor)), nrow=2))$p.value)
chimerism.df$p.adj = p.adjust(chimerism.df$p.value, n = nrow(chimerism.df))

ggplot(chimerism.df, aes(x=donor.freq, y=recipient.freq, color=cluster)) + 
  geom_abline(slope = 1) +
  geom_point() + 
  scale_color_manual(values = Tcell.colors) +
  scale_x_log10('%T cells donor', limits = c(0.001, 0.5), breaks = c(0.001, 0.01, 0.1, 1), labels = c('ND', '1%', '10%', '100%')) + 
  scale_y_log10('%T cells recipient',limits = c(0.001, 0.5), breaks = c(0.001, 0.01, 0.1, 1), labels = c('ND', '1%', '10%', '100%')) + 
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.title = element_text('Arial', size=10, colour = 'black'),
        axis.text = element_text('Arial', size=10, colour = 'black'))
ggsave('./figures/plots/20210610_1026_donor_recipient.svg', width = 2.5, height = 2.5)

# T cells donor versus recipient
donor.cells = cells.ordered[which(chimerism.analysis == 'donor')]
recipient.cells = cells.ordered[which(chimerism.analysis == 'recipient')]
DimPlot(TSB.Tcell, cells.highlight = recipient.cells) + NoLegend()

summary(rowMeans(combined.frequencies[donor.variants, donor.cells]))
summary(rowMeans(combined.frequencies[recipient.variants, recipient.cells]))

# chimerism changes
chimerism.df = c(sample = as.character(), cells.sample = as.numeric(), donor.cells = as.numeric(), recipient.cells = as.numeric())
for (s in c('AML1026_1', 'AML1026_2', 'AML1026_3', 'AML1026_4')) {
  cells.sample = AML.Tcell$cellNames[which(AML.Tcell$Sample.hash == s)]
  donor.cells = length(which(chimerism.analysis[which(cells.ordered %in% cells.sample)] == 'donor'))
  recipient.cells = length(which(chimerism.analysis[which(cells.ordered %in% cells.sample)] == 'recipient'))
  chimerism.df = rbind(chimerism.df, data.frame(sample = s, cells.sample = length(cells.sample), 
                                                donor.cells = donor.cells, recipient.cells = recipient.cells))
}
chimerism.df$chimerism = chimerism.df$donor.cells / (chimerism.df$donor.cells + chimerism.df$recipient.cells)


min(rowMeans(combined.frequencies[donor.variants, donor.cells]))
max(rowMeans(combined.frequencies[donor.variants, donor.cells]))
min(rowMeans(combined.frequencies[recipient.variants, recipient.cells]))
max(rowMeans(combined.frequencies[recipient.variants, recipient.cells]))


