# mtDNA mutation analysis of T cells in AML1012
# used for ASH 2021 abstract (https://doi.org/10.1182/blood-2021-145304)

library(ArchR)
library(gplots)
library(Seurat)
library(ComplexHeatmap)

source('./analysis/Tcell.colors.R')

AML.1012.mito = loadArchRProject('./data/10026/AML.1012.mito/')
AML.Tcell = loadArchRProject('./data/10026/AML.Tcell/')

### 1012
s = 'AML1012'
donor.variants = gtools::mixedsort(c('10463T>C', '16519T>C', '15928G>A', '16126T>C', '16153G>A', '73A>G', '15607A>G', '15452C>A', '7028C>T', '4917A>G',
                                     '13368G>A', '4216T>C', '8697G>A', '11812A>G', '14766C>T', '11719G>A', '709G>A', '8269G>A', '14905G>A', '2706A>G',
                                     '11251A>G', '14233A>G', '1888G>A', '9947G>A', '150C>T', '16294C>T', '16296C>T'))
recipient.variants = gtools::mixedsort(c('16304T>C', '456C>T', '8433T>C', '15833C>T', '4336T>C', '9722T>C', '4011C>T', '93A>G'))

# read mtDNA data
variants = read.table(paste0('./data/10026/mtDNA/20210429_',s,'_high_confidence_variants.csv'), header = T)
combined.frequencies = readRDS(paste0('./data/10026/mtDNA/20210429_',s,'_combined_mutation_frequencies.rds'))

# cells with mtDNA information
cells = AML.Tcell$cellNames[which(grepl(s,AML.Tcell$Sample))]
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

# filter uninformative mtDNA mutations
combined.frequencies = combined.frequencies[-which(rowSums(combined.frequencies) == 0),]
combined.frequencies = combined.frequencies[which(rowMeans(combined.frequencies) > 0.001),]
#combined.frequencies = combined.frequencies[variants$variant[which(variants$n_cells_conf_detected > 50)],]
#combined.frequencies = combined.frequencies[order(rowVars(as.matrix(combined.frequencies)), decreasing = T)[1:30],]

boo = AML.1012.mito$cellNames[which(AML.1012.mito$manual.clusters %in% c('CD4', 'CD8'))]
boo = boo[which(boo %in% colnames(combined.frequencies))]

mtDNA.variants = names(sort(rowSums(combined.frequencies[,boo]), decreasing = T))[1:20]

# order by T cell phenotype
cells.ordered = c()
for (cluster in names(Tcell.colors)) {
  boo = AML.Tcell$cellNames[which(AML.Tcell$manual_clusters == cluster & AML.Tcell$cellNames %in% cells)]
  names(boo) = rep(cluster, length(boo))
  cells.ordered = c(cells.ordered, boo)
}
cells.ordered = cells.ordered[-which(colSums(combined.frequencies[mtDNA.variants, cells.ordered]) == 0)]

# score calculated from donor and recipient variants
chimerism.score = as.numeric(colMeans(combined.frequencies[donor.variants,cells.ordered])) - 
  as.numeric(colMeans(combined.frequencies[recipient.variants,cells.ordered]))

boo=as.data.frame(chimerism.score)
ggplot(data=boo, aes(x=chimerism.score)) + geom_histogram()

chimerism.analysis = ifelse(chimerism.score > 0, 'donor', 'recipient')

source ('./R/20200328_mtscatac_seq.R')
cluster.information = cluster_relevant_mutations(combined.frequencies[mtDNA.variants, cells.ordered], resolution = 1.5, dims = 8, k_param = 30)
cells.ordered.2 = names(cluster.information[[1]])

cluster.df = data.frame(cellbarcode = AML.Tcell$cellNames, manual.cluster = AML.Tcell$manual_clusters)
rownames(cluster.df) = cluster.df$cellbarcode

# plot heatmap
ha.top = ComplexHeatmap::columnAnnotation(manual_cluster = cluster.df[cells.ordered.2, 'manual.cluster'],
                                          #chimerism = chimerism.analysis,
                                          show_legend=F,show_annotation_name = T, 
                                          #annotation_label = c('manual_cluster' = 'phenotype', 'chimerism' = 'origin'),
                                          annotation_name_gp = gpar(fontsize= 8),
                                          col = list(manual_cluster = Tcell.colors,
                                                     chimerism = c('donor' = 'darkblue', 'recipient' = 'darkred')), border=T)

ha.left = ComplexHeatmap::rowAnnotation(chimerism = c(rep('donor', length(donor.variants)), rep('recipient', length(recipient.variants))), 
                                        show_legend =F, show_annotation_name = F, 
                                        col = list(chimerism = c('donor' = 'darkblue', 'recipient' = 'darkred')))

col_fun = circlize::colorRamp2(seq(0.001, 0.1, length=100), colorRampPalette(BuenColors::jdb_palette(name = 'solar_rojos', n=9))(100))

mtDNA.variants = c("1608G>C", "5570T>C", "1970G>C", "779T>C",  "740G>C","3063G>C", "16187C>G","5237G>A",  "709G>A", "16390G>A", "2636G>C", "13759G>C",
                   "2573G>A","3010G>C","6931G>C",  "15731G>C", "1147G>C", "3352G>C",  "11456G>C", "15549T>C")

svglite::svglite('./figure_10026/figures/heatmaps/20210607_AML1012_mtDNA.svg', width = 3, height = 3)
ComplexHeatmap::Heatmap(as.matrix(combined.frequencies[mtDNA.variants,cells.ordered.2]), 
                        show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F,
                        show_column_names = F, show_row_names = T, show_heatmap_legend = F, 
                        row_title = '20 mtDNA mutations', row_title_gp = gpar(fontsize=10, color='black'),
                        column_title = paste0(length(cells.ordered.2), ' T cells'), 
                        column_title_side = 'bottom', column_title_gp = gpar(fontsize=10, color='black'),
                        top_annotation = ha.top, column_split = cluster.information[[2]],
                        #left_annotation = ha.left, 
                        col = col_fun, 
                        row_names_gp = gpar(fontsize = 8), use_raster = T, raster_quality = 10, border=T)
dev.off()

saveRDS(file='./data/10026/objects/20211107_cluster_information_AML1012_Tcell.rds', cluster.information)