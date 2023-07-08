# donor-recipient deconvolution of sample CLL2

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)

col_fun = circlize::colorRamp2(breaks = seq(0,100,100/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))

vafs.1 = as.data.frame(data.table::fread('./data/coevolution/CLL2/AF.csv'))
vafs.1$Barcode = paste0('CLL2#',vafs.1$Barcode)
rownames(vafs.1) = vafs.1$Barcode
vafs.1 = vafs.1[,-c(1,2)]

# explore SNPs
boo = vafs.1[sample(nrow(vafs.1), size = 1000),]
Heatmap(t(boo), show_row_names = T, show_column_names = F, show_row_dend = F, show_column_dend = F, row_names_side = 'left', 
        row_names_gp = gpar(fontsize=8),
        column_km = 2, use_raster = T, raster_quality = 10,
        col = col_fun, border = T, row_split = ifelse(grepl('chrM',colnames(boo)), 'mito','auto'))

# recipient SNPs
variants.CLL3.recipient = c('ASL:chr7:65554385:C/T', 'ASL:chr7:65554352:C/T', 'RELA:chr11:65429479:C/G', 'B3GALT1:chr2:168726231:C/T', 'ESRRG:chr1:216880662:G/A',
                            'GRIK5:chr19:42546856:C/T', 'UPP2:chr2:158974339:G/A', 'TINAG:chr6:54191662:G/A', 'ALDH1A1:chr9:75543910:G/A', 'MME:chr3:154866416:G/T',
                            'PELI2:chr14:56645144:A/T', 'KCNK13:chr14:90650709:G/A', 'CR1:chr1:207791577:G/A', 'chr4:155163826:A/T', 'ABCB5:chr7:20782612:T/C',
                            'ASL:chr7:65554306:A/G', 'U2AF2:chr19:56172500:T/A')
# donor SNPs
variants.CLL3.donor = c('TINAG:chr6:54191686:C/T', 'chr4:155163784:TCA/T')

df = data.frame(variants.CLL3.recipient = rowMeans(vafs.1[,variants.CLL3.recipient]),
                variants.CLL3.donor = rowMeans(vafs.1[,variants.CLL3.donor]))

cells.CLL3.recipient = rownames(df)[which(df$variants.CLL3.recipient > 5 & df$variants.CLL3.donor < 5)]
cells.CLL3.donor = rownames(df)[which(df$variants.CLL3.recipient < 5 & df$variants.CLL3.donor > 5)]

annotation = data.frame(bc = rownames(vafs.1), 
                        annotation = NA)
rownames(annotation) = annotation$bc
annotation[cells.CLL3.recipient, 'annotation'] = 'CLL3.recipient'
annotation[cells.CLL3.donor, 'annotation'] = 'CLL3.donor'
annotation$sample = stringr::str_split_fixed(annotation$annotation, pattern = '\\.', n=2)[,1]
annotation$sample[which(annotation$sample == 'CLL3')] = 'CLL3_2'
annotation$individual = stringr::str_split_fixed(annotation$annotation, pattern = '\\.', n=2)[,2]

write.csv2(annotation, file = './data/coevolution/CLL2/20220506_deconvolution.csv', quote = F)

