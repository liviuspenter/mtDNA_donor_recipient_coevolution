# explore somatic nuclear and mtDNA mutations together

library(ComplexHeatmap)
library(ggplot2)
library(ggExtra)
library(Seurat)

source('./R/filehandling.functions.tapestri.R')

# read data using helper functions
AML1026.1.vafs = read_vafs('./data/coevolution/1026_1/AF.csv', 'AML1026.1')
AML1026.3.vafs = read_vafs('./data/coevolution/1026_3/AF.csv', 'AML1026.3')

AML1026.1.vafs.combined = merge_with_mgatk(AML1026.1.vafs, 
                                           './data/coevolution/vafs/20220211_AML1026.1_mtDNA_vafs.csv', 
                                           'AML1026.1', truncate_barcodes = T)
AML1026.3.vafs.combined = merge_with_mgatk(AML1026.3.vafs, 
                                           './data/coevolution/vafs/20220211_AML1026.3_mtDNA_vafs.csv', 
                                           'AML1026.3', truncate_barcodes = T)

# consider mtDNA and nuclear single nucleotide variants detectable across both samples
AML1026.vafs = rbind(AML1026.1.vafs[,intersect(colnames(AML1026.1.vafs), colnames(AML1026.3.vafs))], 
                     AML1026.3.vafs[,intersect(colnames(AML1026.1.vafs), colnames(AML1026.3.vafs))])
AML1026.vafs.combined = rbind(AML1026.1.vafs.combined[,intersect(colnames(AML1026.1.vafs.combined), colnames(AML1026.3.vafs.combined))], 
                              AML1026.3.vafs.combined[,intersect(colnames(AML1026.1.vafs.combined), colnames(AML1026.3.vafs.combined))])

###
# donor-recipient deconvolution using germline SNPs and maternal mtDNA variants
###
variants.1 = c('NF1:chr17:29559932:C/A', 'TP53:chr17:7579801:G/C', 'TP53:chr17:7578115:T/C', '16304T>C', '4336T>C')
variants.2 = c('SF3B1:chr2:198267770:G/GAA',
               'NPM1:chr5:170837457:A/G', 'KDM6A:chrX:44938563:G/A', 'DNMT3A:chr2:25463483:G/A',
               'BRAF:chr7:140449071:C/G', 'FLT3:chr13:28592546:T/C', 'NF1:chr17:29483195:G/C', '16294C>T', '16296C>T')
#RHP.mutations = c('ASXL1'='ASXL1:chr20:31022441:A/G', 'NRAS' = 'NRAS:chr1:115258745:C/G', 'SF3B1' = 'SF3B1:chr2:198266512:C/T')
RHP.mutations = c('NRAS' = 'NRAS:chr1:115258745:C/G', 'SF3B1' = 'SF3B1:chr2:198266512:C/T')
DP1 = as.data.frame(data.table::fread('./data/coevolution/1026_1/DP.csv'))
DP1$Barcode = paste0('AML1026.1#', DP1$Barcode)
DP1[DP1 == 0] = NA
DP3 = as.data.frame(data.table::fread('./data/coevolution/1026_3/DP.csv'))
DP3$Barcode = paste0('AML1026.3#', DP3$Barcode)
DP3[DP3 == 0] = NA
DP = rbind(DP1[c('Sample', 'Barcode', colnames(AML1026.vafs))], DP3[c('Sample', 'Barcode', colnames(AML1026.vafs))])
cells.plot = DP$Barcode[which(complete.cases(DP[,as.character(RHP.mutations)]))]

VARIANT_CUTOFF = 10
variant.df = data.frame(variants.1 = rowMeans(AML1026.vafs.combined[,variants.1]),
                        variants.2 = rowMeans(AML1026.vafs.combined[,variants.2]),
                        barcode = rownames(AML1026.vafs))
variant.df$annotation = 'none'
variant.df$annotation[which(variant.df$variants.1 > VARIANT_CUTOFF & variant.df$variants.2 > VARIANT_CUTOFF)] = 'doublet'
variant.df$annotation[which(variant.df$variants.1 > VARIANT_CUTOFF & variant.df$variants.2 < VARIANT_CUTOFF)] = 'variant1'
variant.df$annotation[which(variant.df$variants.1 < VARIANT_CUTOFF & variant.df$variants.2 > VARIANT_CUTOFF)] = 'variant2'

variant.df = variant.df[cells.plot,]

ha = HeatmapAnnotation(variant = c(variant.df$annotation[which(variant.df$annotation == 'variant1')], 
                                   variant.df$annotation[which(variant.df$annotation == 'variant2')]),
                       col = list(variant = c('variant1' = 'orange', 'variant2' = 'purple')),
                       annotation_legend_param = list(variant = list(title = 'Individual')), border = T)
col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = BuenColors::jdb_palette(name = 'solar_rojos', n=9))

# visualize donor-recipient deconvolution using germline SNPs and maternal mtDNA variants
svglite::svglite('./figure_coevolution/AML/figures/20220212_AML1026_variants_heatmap.svg', width = 5, height = 4)
Heatmap(t(AML1026.vafs.combined[c(variant.df$barcode[which(variant.df$annotation == 'variant1')], 
                 variant.df$barcode[which(variant.df$annotation == 'variant2')]),c(variants.1, variants.2)]), 
        top_annotation = ha, raster_quality = 10, use_raster = T, border = T,
        cluster_rows = T, cluster_columns = F, show_column_names = F, show_row_names = T, 
        show_column_dend = F, show_row_dend = F, row_names_gp = gpar(fontsize=8))
dev.off()

###
# track somatic nuclear (NRAS, SF3B1) and mtDNA mutation (11736T>C) across different recipient-derived cell types 
###
cluster.colors = c('HSC' = RColorBrewer::brewer.pal(name = 'YlOrRd', n = 5)[5],
                   'Progenitor' = RColorBrewer::brewer.pal(name = 'YlOrRd', n = 5)[4],
                   'Mono' = RColorBrewer::brewer.pal(name = 'YlOrRd', n = 5)[3],
                   'Erythroid' = RColorBrewer::brewer.pal(name = 'YlOrRd', n = 5)[2],
                   'CD4' = 'lightblue', 
                   'CD8' = 'darkblue', 
                   'Plasma' = 'darkgreen')

so.13 = readRDS('./data/coevolution/objects/20220112_AML1026.rds')
so.13 = RenameCells(so.13, new.names = paste0('AML', colnames(so.13)))
variant.df$manual.cluster = so.13$manual.cluster[rownames(variant.df)]

vafs = AML1026.vafs.combined[variant.df$barcode[which(variant.df$annotation == 'variant1')],
                             c(RHP.mutations, colnames(AML1026.vafs.combined)[which(grepl('>',colnames(AML1026.vafs.combined)))])]
vafs[, which(!colnames(vafs) %in% RHP.mutations)] = 10*vafs[, which(!colnames(vafs) %in% RHP.mutations)] 

col_fun = circlize::colorRamp2(breaks = seq(0,100,100/8), colors = BuenColors::jdb_palette(name = 'solar_rojos', n=9))

svglite::svglite('./figure_coevolution/AML/figures/20220212_AML1026_mutations_heatmap.svg', width = 5, height = 1.5)
ha = HeatmapAnnotation(celltype = variant.df[rownames(vafs), 'manual.cluster'], col = list(celltype = cluster.colors),
                       border = T, simple_anno_size = unit(5, 'pt'), annotation_name_gp = gpar(fontsize=8), gp = gpar(fontsize=8),
                       annotation_legend_param = )
Heatmap(t(vafs[,unique(c(rev(RHP.mutations), '11736T>C'))]), top_annotation = ha, raster_quality = 10, use_raster = T, border = T,
        cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T, 
        show_column_dend = F, show_row_dend = F, row_names_gp = gpar(fontsize=8), column_split = variant.df[rownames(vafs), 'manual.cluster'],
        col = col_fun, column_title = paste0(nrow(vafs), ' cells'), column_title_side = 'bottom', column_title_gp = gpar(fontsize=8))
dev.off()

# visualize donor-recipient deconvolution on UMAP
p=DimPlot(so.13, group.by = 'annotation', cols = c('none' = 'grey', 'doublet' = 'black', 
                                                 'variant1' = 'purple', 'variant2' = 'orange')) + 
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave('./figure_coevolution/AML/figures/20230428_donor_recipient_UMAP.png', width = 4, height = 4, dpi = 600, plot = p)
