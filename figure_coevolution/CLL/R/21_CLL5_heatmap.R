# create heatmap of CLL5 (called CLL3)

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)

`%ni%` <- Negate(`%in%`)

col_fun = circlize::colorRamp2(breaks = c(-1,seq(0,100,100/8)), colors = c('grey', BuenColors::jdb_palette(name = 'solar_rojos')))

# call mtDNA mutations using mgatk 
CLL2.mito = ReadMGATK('./data/coevolution/CLL2.hg38.mgatk/')
colnames(CLL2.mito$counts) = gsub(paste0('CLL2#', colnames(CLL2.mito$counts)), pattern = '-1', replacement = '')
rownames(CLL2.mito$depth) = gsub(paste0('CLL2#', rownames(CLL2.mito$depth)), pattern = '-1', replacement = '')
CLL3.mito = ReadMGATK('./data/coevolution/CLL3.mgatk/')
colnames(CLL3.mito$counts) = gsub(paste0('CLL3#', colnames(CLL3.mito$counts)), pattern = '-1', replacement = '')
rownames(CLL3.mito$depth) = gsub(paste0('CLL3#', rownames(CLL3.mito$depth)), pattern = '-1', replacement = '')
CLL.mito = CLL2.mito
CLL.mito$depth = rbind(CLL.mito$depth, CLL3.mito$depth)
CLL.mito$counts = cbind(CLL.mito$counts, CLL3.mito$counts)

CLL.so <- CreateAssayObject(counts = CLL.mito$counts)
CLL.mgatk = CreateSeuratObject(counts = CLL.so, project = 'CLL.mgatk', assay = 'mito')
rm(CLL.so)
CLL.mgatk = AddMetaData(CLL.mgatk, metadata = CLL.mito$depth, col.name = 'mtDNA_depth')

CLL.mgatk <- AlleleFreq(object = CLL.mgatk, variants = c('2332C>T', '5979G>A', '12685T>C', '2776G>A', '4412G>A', 
                                                         '5458T>C', '11886T>C', '13676A>G', '1489G>A', '11700T>C',
                                                         '1263G>A', '279T>C'),
                        assay = "mito")
CLL.mito.vaf = as.data.frame(t(as.data.frame(GetAssayData(CLL.mgatk[["alleles"]]))))

so.12 = readRDS('./data/coevolution/objects/20220504_CLL34_filtered.rds')

vafs.1 = as.data.frame(data.table::fread('./data/coevolution/CLL2.whitelist/AF.csv'))
vafs.1$Barcode = gsub(paste0('CLL2#',vafs.1$Barcode), pattern = '-1', replacement = '')
rownames(vafs.1) = vafs.1$Barcode
vafs.1 = vafs.1[,-c(1,2)]
vafs.2 = as.data.frame(data.table::fread('./data/coevolution/CLL3.whitelist/AF.csv'))
vafs.2$Barcode = gsub(paste0('CLL3#',vafs.2$Barcode), pattern = '-1', replacement = '')
rownames(vafs.2) = vafs.2$Barcode
vafs.2 = vafs.2[,-c(1,2)]

depth.1 = as.data.frame(data.table::fread('./data/coevolution/CLL2.whitelist/DP.csv'))
depth.1$Barcode = gsub(paste0('CLL2#',depth.1$Barcode), pattern = '-1', replacement = '')
rownames(depth.1) = depth.1$Barcode
depth.1 = depth.1[,-c(1,2)]
depth.2 = as.data.frame(data.table::fread('./data/coevolution/CLL3.whitelist/DP.csv'))
depth.2$Barcode = gsub(paste0('CLL3#',depth.2$Barcode), pattern = '-1', replacement = '')
rownames(depth.2) = depth.2$Barcode
depth.2 = depth.2[,-c(1,2)]

vafs = dplyr::bind_rows(vafs.1, vafs.2)
depth = dplyr::bind_rows(depth.1, depth.2)

annotation2 = read.csv2('./data/coevolution/CLL2/20220506_deconvolution.csv', row.names = 1)
annotation3 = read.csv2('./data/coevolution/CLL3/20220506_deconvolution.csv', row.names = 1)

annotation2$bc = gsub(annotation2$bc, pattern = '-1', replacement = '')
annotation3$bc = gsub(annotation3$bc, pattern = '-1', replacement = '')


cells.donor.2 = gsub(annotation2$bc[which(annotation2$annotation == 'CLL3.donor' & annotation2$bc %in% rownames(vafs))],
                     pattern = '-1', replacement = '')
cells.recipient.2 = gsub(annotation2$bc[which(annotation2$annotation == 'CLL3.recipient' & annotation2$bc %in% rownames(vafs))],
                         pattern = '-1', replacement = '')
cells.1 = annotation3$bc[which(annotation3$annotation == 'CLL3.recipient' & annotation3$bc %in% rownames(vafs))]

depth = depth[c(cells.1, cells.donor.2, cells.recipient.2),]
vafs = vafs[c(cells.1, cells.donor.2, cells.recipient.2),]

vafs.auto = vafs[,which(!grepl('chrM', colnames(vafs)))]
vafs.auto[vafs.auto < 10] = 0
depth = depth[rownames(vafs.auto),]
depth = depth[, colnames(vafs.auto)]
vafs.auto[is.na(depth) | depth == 0] = -1

variant.df = data.frame(vmr = -log10(apply(vafs.auto, 2, FUN = function(x){mean(x) / var(x)})),
                        vaf = apply(vafs.auto, 2, FUN = function(x){median(x[which(x != 0)])}),
                        freq = apply(vafs.auto, 2, FUN = function(x){length(which(x != 0))}))

vafs.auto = vafs.auto[intersect(rownames(vafs.auto), rownames(CLL.mito.vaf)),]

vafs.mito = 100*CLL.mito.vaf[rownames(vafs.auto),]

vafs = cbind(vafs.auto, vafs.mito)

relevant.mutations = unique(c(rownames(variant.df)[which(variant.df$vmr > 0 & variant.df$vaf > 25)],
                              colnames(vafs.mito)))

ha = columnAnnotation(cluster = so.12$manual.cluster[c(cells.1, cells.recipient.2, cells.donor.2)],
                      col = list('cluster' = c('CLL' = 'orange', 'CD4 T cell' = 'lightblue', 'CD8 T cell' = 'darkblue', 'NK' = 'purple','Myeloid' = 'darkgreen',
                                               'Cluster 6' = 'black')), 
                      simple_anno_size = unit(5, 'pt'), border = T)

svglite::svglite('./figure_coevolution/CLL/figures/CLL5/heatmaps/20230427_CLL5_raw.svg', width = 15, height = 15)
Heatmap(t(vafs[c(cells.1, cells.recipient.2, cells.donor.2),relevant.mutations]), show_row_names = T, show_column_names = F, 
        show_row_dend = F, show_column_dend = F, row_names_side = 'left', cluster_columns = F, cluster_rows = F,
        row_names_gp = gpar(fontsize=8),
        use_raster = T, raster_quality = 10, column_split = c(rep('pre-FCR', length(cells.1)), 
                                                              rep('post-RIC', length(cells.recipient.2)), 
                                                              rep('donor', length(cells.donor.2))),
        col = col_fun, border = T, row_split = factor(ifelse(grepl('chrM',relevant.mutations), 'mito','auto'), levels = c('auto', 'mito')),
        top_annotation = ha)
dev.off()


### make curated heatmap with custom ordering
source('./R/20200328_mtscatac_seq.R')

vafs = dplyr::bind_rows(vafs.1, vafs.2)
depth = dplyr::bind_rows(depth.1, depth.2)

cells.2.donor = annotation2$bc[which(annotation2$annotation == 'CLL3.donor' & annotation2$bc %in% rownames(vafs))]
cells.2.recipient = annotation2$bc[which(annotation2$annotation == 'CLL3.recipient' & annotation2$bc %in% rownames(vafs))]
cells.1.CLL = annotation3$bc[which(annotation3$annotation == 'CLL3.recipient' & annotation3$bc %in% rownames(vafs) & 
                                     annotation3$bc %in% colnames(so.12)[which(so.12$manual.cluster == 'CLL')])]
cells.1.immune = annotation3$bc[which(annotation3$annotation == 'CLL3.recipient' & annotation3$bc %in% rownames(vafs) & 
                                        annotation3$bc %in% colnames(so.12)[which(so.12$manual.cluster %ni% c('CLL', 'Cluster 6'))])]

CLL.auto = c('B3GALT1:chr2:168726231:C/T', 'GRIK5:chr19:42546856:C/T',
             'TINAG:chr6:54191662:G/A', 'UPP2:chr2:158974339:G/A', 'PELI2:chr14:56645144:A/T', 'ALDH1A1:chr9:75543910:G/A', 'RELA:chr11:65429479:C/G',
             'U2AF2:chr19:56172500:T/A', 'CR1:chr1:207791577:G/A', 'chr4:155163826:A/T', 'ABCB5:chr7:20782612:T/C', 'ASL:chr7:65554306:A/G', 
             'MME:chr3:154866416:G/T', 'KCNK13:chr14:90650709:G/A')
CLL.mito = c('2332C>T', '5979G>A')

depth = depth[c(cells.1.CLL, cells.1.immune, cells.2.recipient, cells.2.donor),]
vafs = vafs[c(cells.1.CLL, cells.1.immune, cells.2.recipient, cells.2.donor),]

vafs.auto = vafs[,c(CLL.auto)]
vafs.auto[vafs.auto < 10] = 0
depth = depth[rownames(vafs.auto),]
depth = depth[, colnames(vafs.auto)]

vafs.auto[is.na(depth) | depth == 0] = NA
vafs.auto = vafs.auto[complete.cases(vafs.auto),]

vafs.auto = vafs.auto[intersect(rownames(vafs.auto), rownames(CLL.mito.vaf)),]

# plot from 0 to 10%
vafs.mito = 1000*CLL.mito.vaf[rownames(vafs.auto),]

vafs = cbind(vafs.auto, vafs.mito)

cells.2.donor = annotation2$bc[which(annotation2$annotation == 'CLL3.donor' & annotation2$bc %in% rownames(vafs))]
cells.2.recipient = annotation2$bc[which(annotation2$annotation == 'CLL3.recipient' & annotation2$bc %in% rownames(vafs))]
cells.1.CLL = annotation3$bc[which(annotation3$annotation == 'CLL3.recipient' & annotation3$bc %in% rownames(vafs) & 
                                     annotation3$bc %in% colnames(so.12)[which(so.12$manual.cluster == 'CLL')])]
cells.1.immune = annotation3$bc[which(annotation3$annotation == 'CLL3.recipient' & annotation3$bc %in% rownames(vafs) & 
                                        annotation3$bc %in% colnames(so.12)[which(so.12$manual.cluster %ni% c('CLL', 'Cluster 6'))])]

cluster.information.1 = cluster_relevant_mutations(t(vafs[cells.1.CLL, c(CLL.auto, CLL.mito)]), dims = 2, k_param = 5)
cells.1.CLL = colnames(cluster.information.1[[1]])
#cells.1.CLL = cells.1.CLL[order(vafs[cells.1.CLL,'ASL:chr7:65554385:C/T'], decreasing = T)]
cluster.information.2 = cluster_relevant_mutations(t(vafs[cells.2.recipient, c(CLL.auto, CLL.mito)]), dims = 2, k_param = 5)
cells.2.recipient = colnames(cluster.information.2[[1]])
#cells.2.recipient = cells.2.recipient[order(vafs[cells.2.recipient,'ASL:chr7:65554385:C/T'], decreasing = T)]
cells.1.immune = names(sort(so.12$manual.cluster[cells.1.immune]))

ha = columnAnnotation(cluster = so.12$manual.cluster[c(cells.1.CLL, cells.1.immune, cells.2.recipient, cells.2.donor)],
                      col = list('cluster' = c('CLL' = 'orange', 'CD4 T cell' = 'lightblue', 'CD8 T cell' = 'darkblue', 'NK' = 'purple','Myeloid' = 'darkgreen',
                                               'Cluster 6' = 'black')), 
                      simple_anno_size = unit(5, 'pt'), border = T)

svglite::svglite('./figure_coevolution/CLL/figures/CLL5/heatmaps/20230427_CLL5_curated.svg', width = 7, height = 2.5)
Heatmap(t(vafs[c(cells.1.CLL, cells.1.immune, cells.2.recipient, cells.2.donor),
               c(CLL.auto, CLL.mito)]), 
        show_row_names = T, show_column_names = F, 
        show_row_dend = F, show_column_dend = F, row_names_side = 'left', cluster_columns = F,
        row_names_gp = gpar(fontsize=8), cluster_rows = F, 
        use_raster = T, raster_quality = 10, 
        row_split = factor(c(rep('CLL auto', length(CLL.auto)),
                             rep('CLL mito', length(CLL.mito))), 
                           levels = c('CLL auto', 'CLL mito')),
        column_split = factor(c(rep('CLL pre-FCR', length(cells.1.CLL)),
                                rep('Immune cells pre-FCR', length(cells.1.immune)),
                                rep('CLL post-RIC', length(cells.2.recipient)), 
                                rep('Donor', length(cells.2.donor))), 
                              levels = c( 'CLL pre-FCR', 'CLL post-RIC', 'Immune cells pre-FCR','Donor')),
        column_title_gp = gpar(fontsize=8), 
        row_title_gp = gpar(fontsize=8),
        col = col_fun, border = T, 
        top_annotation = ha)
dev.off()


### plot individual mutations

ggplot() + 
  geom_point(data=vafs[cells.1.CLL,], aes(x=`MME:chr3:154866416:G/T`, y=`KCNK13:chr14:90650709:G/A`), color='blue', size=0.5) + 
  scale_x_continuous('% MMEK525N', limits = c(0,100)) + 
  scale_y_continuous('% KCKN13V197I', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_pre_FCR_MME_KCKN13.svg', width = 1.1, height = 1.1)

ggplot() + 
  geom_point(data=vafs[cells.2.recipient,], aes(x=`MME:chr3:154866416:G/T`, y=`KCNK13:chr14:90650709:G/A`), color='firebrick', size=0.5) + 
  scale_x_continuous('% MMEK525N', limits = c(0,100)) + 
  scale_y_continuous('% KCKN13V197I', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_post_HSCT_MME_KCKN13.svg', width = 1.1, height = 1.1)


ggplot() + 
  geom_point(data=vafs[cells.1.CLL,], aes(x=`2332C>T`/10, y=`5979G>A`/10), color='blue', size=0.5) + 
  scale_x_continuous('% 2332C>T', limits = c(0,100)) + 
  scale_y_continuous('% 5979G>A', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_pre_FCR_2332C>T_5979G>A.svg', width = 1.1, height = 1.1)

ggplot() + 
  geom_point(data=vafs[cells.2.recipient,], aes(x=`2332C>T`/10, y=`5979G>A`/10), color='firebrick', size=0.5) + 
  scale_x_continuous('% 2332C>T', limits = c(0,100)) + 
  scale_y_continuous('% 5979G>A', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_post_HSCT_2332C>T_5979G>A.svg', width = 1.1, height = 1.1)



ggplot() + 
  geom_point(data=vafs[cells.1.CLL,], aes(x=`MME:chr3:154866416:G/T`, y=`5979G>A`/10), color='blue', size=0.5) + 
  scale_x_continuous('% MMEK525N', limits = c(0,100)) + 
  scale_y_continuous('% 5979G>A', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_pre_FCR_MME_5979G>A.svg', width = 1.1, height = 1.1)

ggplot() + 
  geom_point(data=vafs[cells.2.recipient,], aes(x=`MME:chr3:154866416:G/T`, y=`5979G>A`/10), color='firebrick', size=0.5) + 
  scale_x_continuous('% MMEK525N', limits = c(0,100)) + 
  scale_y_continuous('% 5979G>A', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_post_HSCT_MME_5979G>A.svg', width = 1.1, height = 1.1)


ggplot() + 
  geom_point(data=vafs[cells.1.CLL,], aes(x=`MME:chr3:154866416:G/T`, y=`2332C>T`/10), color='blue', size=0.5) + 
  scale_x_continuous('% MMEK525N', limits = c(0,100)) + 
  scale_y_continuous('% 2332C>T', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_pre_FCR_MME_2332C>T.svg', width = 1.1, height = 1.1)

ggplot() + 
  geom_point(data=vafs[cells.2.recipient,], aes(x=`MME:chr3:154866416:G/T`, y=`2332C>T`/10), color='firebrick', size=0.5) + 
  scale_x_continuous('% MMEK525N', limits = c(0,100)) + 
  scale_y_continuous('% 2332C>T', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_post_HSCT_MME_2332C>T.svg', width = 1.1, height = 1.1)



ggplot() + 
  geom_point(data=vafs[cells.1.CLL,], aes(x=`ASL:chr7:65554306:A/G`, y=`2332C>T`/10), color='blue', size=0.5) + 
  scale_x_continuous('% ASLY321C', limits = c(0,100)) + 
  scale_y_continuous('% 2332C>T', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_pre_FCR_ASL_2332C>T.svg', width = 1.1, height = 1.1)

ggplot() + 
  geom_point(data=vafs[cells.2.recipient,], aes(x=`ASL:chr7:65554306:A/G`, y=`2332C>T`/10), color='firebrick', size=0.5) + 
  scale_x_continuous('% ASLY321C', limits = c(0,100)) + 
  scale_y_continuous('% 2332C>T', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_post_HSCT_ASL_2332C>T.svg', width = 1.1, height = 1.1)

ggplot() + 
  geom_point(data=vafs[cells.1.CLL,], aes(x=`ASL:chr7:65554306:A/G`, y=`5979G>A`/10), color='blue', size=0.5) + 
  scale_x_continuous('% ASLY321C', limits = c(0,100)) + 
  scale_y_continuous('% 5979G>A', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_pre_FCR_ASL_5979G>A.svg', width = 1.1, height = 1.1)

ggplot() + 
  geom_point(data=vafs[cells.2.recipient,], aes(x=`ASL:chr7:65554306:A/G`, y=`5979G>A`/10), color='firebrick', size=0.5) + 
  scale_x_continuous('% ASLY321C', limits = c(0,100)) + 
  scale_y_continuous('% 5979G>A', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_post_HSCT_ASL_5979G>A.svg', width = 1.1, height = 1.1)




ggplot() + 
  geom_point(data=vafs[cells.1.CLL,], aes(x=`ASL:chr7:65554306:A/G`, y=`KCNK13:chr14:90650709:G/A`), color='blue', size=0.5) + 
  scale_x_continuous('% ASLY321C', limits = c(0,100)) + 
  scale_y_continuous('% KCKN13V197I', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_pre_FCR_ASL_KCNK13.svg', width = 1.1, height = 1.1)

ggplot() + 
  geom_point(data=vafs[cells.2.recipient,], aes(x=`ASL:chr7:65554306:A/G`, y=`KCNK13:chr14:90650709:G/A`), color='firebrick', size=0.5) + 
  scale_x_continuous('% ASLY321C', limits = c(0,100)) + 
  scale_y_continuous('% KCKN13V197I', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_post_HSCT_ASL_KCNK13.svg', width = 1.1, height = 1.1)

ggplot() + 
  geom_point(data=vafs[cells.1.CLL,], aes(x=`ASL:chr7:65554306:A/G`, y=`MME:chr3:154866416:G/T`), color='blue', size=0.5) + 
  scale_x_continuous('% ASLY321C', limits = c(0,100)) + 
  scale_y_continuous('% MMEK525N', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_pre_FCR_ASL_MME.svg', width = 1.1, height = 1.1)

ggplot() + 
  geom_point(data=vafs[cells.2.recipient,], aes(x=`ASL:chr7:65554306:A/G`, y=`MME:chr3:154866416:G/T`), color='firebrick', size=0.5) + 
  scale_x_continuous('% ASLY321C', limits = c(0,100)) + 
  scale_y_continuous('% MMEK525N', limits = c(0,100)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./figure_coevolution/CLL/figures/CLL5/plots/20230427_CLL5_post_HSCT_ASL_MME.svg', width = 1.1, height = 1.1)
