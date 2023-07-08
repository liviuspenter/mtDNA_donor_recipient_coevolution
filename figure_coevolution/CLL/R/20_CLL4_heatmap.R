# create heatmap of CLL4 (called CLL1)

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)

`%ni%` <- Negate(`%in%`)

col_fun = circlize::colorRamp2(breaks = c(-1,seq(0,100,100/8)), colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))

# call mtDNA mutations using mgatk from samples CLL3 and CLL4
CLL3.mito = ReadMGATK('./data/coevolution/CLL3.mgatk/')
colnames(CLL3.mito$counts) = gsub(paste0('CLL3#', colnames(CLL3.mito$counts)), pattern = '-1', replacement = '')
rownames(CLL3.mito$depth) = gsub(paste0('CLL3#', rownames(CLL3.mito$depth)), pattern = '-1', replacement = '')
CLL4.mito = ReadMGATK('./data/coevolution/CLL4.mgatk/')
colnames(CLL4.mito$counts) = gsub(paste0('CLL4#', colnames(CLL4.mito$counts)), pattern = '-1', replacement = '')
rownames(CLL4.mito$depth) = gsub(paste0('CLL4#', rownames(CLL4.mito$depth)), pattern = '-1', replacement = '')
CLL.mito = CLL3.mito
CLL.mito$depth = rbind(CLL.mito$depth, CLL4.mito$depth)
CLL.mito$counts = cbind(CLL.mito$counts, CLL4.mito$counts)

CLL.so <- CreateAssayObject(counts = CLL.mito$counts)
CLL.mgatk = CreateSeuratObject(counts = CLL.so, project = 'CLL.mgatk', assay = 'mito')
rm(CLL.so)
CLL.mgatk = AddMetaData(CLL.mgatk, metadata = CLL.mito$depth, col.name = 'mtDNA_depth')

CLL.mgatk <- AlleleFreq(object = CLL.mgatk, variants = c('3538G>A', '16247A>G', '6426G>A', '16290C>T',
                                                         '13916G>A', '4344T>C', '5650G>A', '786G>A', '1918G>A'),
                         assay = "mito")
CLL.mito.vaf = as.data.frame(t(as.data.frame(GetAssayData(CLL.mgatk[["alleles"]]))))

so.12 = readRDS('./data/coevolution/objects/20220504_CLL34_filtered.rds')

annotation3 = read.csv2('./data/coevolution/CLL3/20220506_deconvolution.csv', row.names = 1)
annotation4 = read.csv2('./data/coevolution/CLL4/20220506_deconvolution.csv', row.names = 1)

cells.donor.2 = annotation3$bc[which(annotation3$annotation == 'CLL1.donor')]
cells.recipient.2 = annotation3$bc[which(annotation3$annotation == 'CLL1.recipient')]
cells.1 = annotation4$bc[which(annotation4$annotation == 'CLL1.recipient')]

vafs.1 = as.data.frame(data.table::fread('./data/coevolution/CLL3.whitelist/AF.csv'))
vafs.1$Barcode = paste0('CLL3#',vafs.1$Barcode)
rownames(vafs.1) = vafs.1$Barcode
vafs.1 = vafs.1[,-c(1,2)]
vafs.2 = as.data.frame(data.table::fread('./data/coevolution/CLL4.whitelist/AF.csv'))
vafs.2$Barcode = paste0('CLL4#',vafs.2$Barcode)
rownames(vafs.2) = vafs.2$Barcode
vafs.2 = vafs.2[,-c(1,2)]

depth.1 = as.data.frame(data.table::fread('./data/coevolution/CLL3.whitelist/DP.csv'))
depth.1$Barcode = paste0('CLL3#',depth.1$Barcode)
rownames(depth.1) = depth.1$Barcode
depth.1 = depth.1[,-c(1,2)]
depth.2 = as.data.frame(data.table::fread('./data/coevolution/CLL4.whitelist/DP.csv'))
depth.2$Barcode = paste0('CLL4#',depth.2$Barcode)
rownames(depth.2) = depth.2$Barcode
depth.2 = depth.2[,-c(1,2)]

depth = dplyr::bind_rows(depth.1, depth.2)
depth = depth[c(cells.1, cells.donor.2, cells.recipient.2),]

vafs = dplyr::bind_rows(vafs.1, vafs.2)
vafs = vafs[c(cells.1, cells.donor.2, cells.recipient.2),]
vafs.auto = vafs[,which(!grepl('chrM', colnames(vafs)))]
vafs.auto[vafs.auto < 10] = 0
depth = depth[rownames(vafs.auto),]
depth = depth[, colnames(vafs.auto)]
vafs.auto[is.na(depth) | depth == 0] = -1

variant.df = data.frame(vmr = -log10(apply(vafs.auto, 2, FUN = function(x){mean(x) / var(x)})),
                        vaf = apply(vafs.auto, 2, FUN = function(x){median(x[which(x != 0)])}),
                        freq = apply(vafs.auto, 2, FUN = function(x){length(which(x != 0))}))

vafs.mito = 100*CLL.mito.vaf[rownames(vafs.auto),]

vafs = cbind(vafs.auto, vafs.mito)

relevant.mutations = unique(c(rownames(variant.df)[which(variant.df$vmr > 0 & variant.df$vaf > 25)],
                              colnames(vafs.mito)))

ha = columnAnnotation(cluster = so.12$manual.cluster[c(cells.1, cells.recipient.2, cells.donor.2)],
                      col = list('cluster' = c('CLL' = 'orange', 'CD4 T cell' = 'lightblue', 'CD8 T cell' = 'darkblue', 'NK' = 'purple','Myeloid' = 'darkgreen',
                                               'Cluster 6' = 'orange')), 
                      simple_anno_size = unit(5, 'pt'), border = T)

svglite::svglite('./figure_coevolution/CLL/figures/CLL4/heatmaps/20230426_CLL4_raw.svg', width = 15, height = 15)
Heatmap(t(vafs[c(cells.1, cells.recipient.2, cells.donor.2),relevant.mutations]), show_row_names = T, show_column_names = F, 
        show_row_dend = F, show_column_dend = F, row_names_side = 'left', cluster_columns = T,
        row_names_gp = gpar(fontsize=8),
        use_raster = T, raster_quality = 10, column_split = c(rep('pre-FCR', length(cells.1)), 
                                                              rep('post-RIC', length(cells.recipient.2)), 
                                                              rep('donor', length(cells.donor.2))),
        col = col_fun, border = T, row_split = factor(ifelse(grepl('chrM',relevant.mutations), 'mito','auto'), levels = c('auto', 'mito')),
        top_annotation = ha)
dev.off()

### make curated heatmap with custom ordering
source('./R/20200328_mtscatac_seq.R')
donor.auto = c('TINAG:chr6:54191736:C/T')
donor.mito = c('1918G>A','5650G>A', '786G>A')
recipient.auto = c('BTK:chrX:100611285:T/C')
CLL.common.auto = c('BCAM:chr19:45317870:G/A', 'PPM1E:chr17:57043202:T/C', 'SF3B1:chr2:198267371:G/T', 'SH3GL1:chr19:4363776:T/A', 'WNK1:chr12:966369:G/T')
CLL.common.mito = c('3538G>A')
CLL.clone1.auto = c('OPRK1:chr8:54163591:A/T', 'IGLL5:chr22:23235972:G/A', 'UBR4:chr1:19439249:C/T', 'ICE1:chr5:5464776:C/T', 'CA3:chr8:86360275:C/T', 
                    'NAA40:chr11:63721451:A/T', 'LAMA5:chr20:60912714:C/A')
CLL.clone1.mito = c('4344T>C', '13916G>A', '16247A>G')
CLL.clone2.auto = c('KCNA6:chr12:4920442:C/T', 'SH3RF1:chr4:170043260:C/T', 'DPCD:chr10:103354473:C/A', 'ZNF215:chr11:6953628:A/C', 'MGA:chr15:42019401:C/T',
                    'ALPI:chr2:233323457:C/T', 'DOCK1:chr10:128841379:G/A', 'GFOD2:chr16:67709195:C/T', 'MAP7:chr6:136682295:G/A')
CLL.clone2.mito = c('6426G>A', '16290C>T')

cells.2.donor = annotation3$bc[which(annotation3$annotation == 'CLL1.donor' & annotation3$bc %in% rownames(vafs) &
                                             annotation3$bc %in% colnames(so.12)[which(so.12$manual.cluster %ni% c('CLL', 'Cluster 6'))])]
cells.2.recipient = annotation3$bc[which(annotation3$annotation == 'CLL1.recipient' & annotation3$bc %in% rownames(vafs) &
                                                 annotation3$bc %in% colnames(so.12)[which(so.12$manual.cluster %in% c('CLL', 'Cluster 6'))])]
cells.1.CLL = annotation4$bc[which(annotation4$annotation == 'CLL1.recipient' & annotation4$bc %in% rownames(vafs) & 
                                           annotation4$bc %in% colnames(so.12)[which(so.12$manual.cluster == 'CLL')])]
cells.1.immune = annotation4$bc[which(annotation4$annotation == 'CLL1.recipient' & annotation4$bc %in% rownames(vafs) & 
                                              annotation4$bc %in% colnames(so.12)[which(so.12$manual.cluster %ni% c('CLL', 'Cluster 6'))])]
cells.2.immune = annotation3$bc[which(annotation3$annotation == 'CLL1.recipient' & annotation3$bc %in% rownames(vafs) & 
                                              annotation3$bc %in% colnames(so.12)[which(so.12$manual.cluster %ni% c('CLL', 'Cluster 6'))])]

cluster.information.1 = cluster_relevant_mutations(t(vafs[cells.1.CLL, c(CLL.clone1.auto, CLL.clone1.mito)]), dims = 5, k_param = 5)
cells.1.CLL = colnames(cluster.information.1[[1]])
cluster.information.2 = cluster_relevant_mutations(t(vafs[cells.2.recipient, c(CLL.clone2.auto, CLL.clone2.mito)]), dims = 5, k_param = 5)
cells.2.recipient = colnames(cluster.information.2[[1]])
cells.2.donor = cells.2.donor[order(vafs[cells.2.donor,donor.mito], decreasing = T)]
cells.2.donor = names(sort(so.12$manual.cluster[cells.2.donor]))
cells.1.immune = names(sort(so.12$manual.cluster[cells.1.immune]))
cells.2.immune = names(sort(so.12$manual.cluster[cells.2.immune]))

ha = columnAnnotation(cluster = so.12$manual.cluster[c(cells.1.CLL, cells.1.immune, cells.2.recipient, cells.2.donor)],
                      col = list('cluster' = c('CLL' = 'orange', 'CD4 T cell' = 'lightblue', 'CD8 T cell' = 'darkblue', 'NK' = 'purple','Myeloid' = 'darkgreen',
                                               'Cluster 6' = 'orange')), 
                      simple_anno_size = unit(5, 'pt'), border = T)

svglite::svglite('./figure_coevolution/CLL/figures/CLL4/heatmaps/20230426_CLL4_curated.svg', width = 6, height = 3.2)
Heatmap(t(vafs[c(cells.1.CLL, cells.1.immune, cells.2.recipient, cells.2.donor),
               c(CLL.common.auto, CLL.common.mito, CLL.clone1.auto, CLL.clone1.mito, CLL.clone2.auto, CLL.clone2.mito, donor.mito)]), 
        show_row_names = T, show_column_names = F, 
        show_row_dend = F, show_column_dend = F, row_names_side = 'left', cluster_columns = F,
        row_names_gp = gpar(fontsize=8), cluster_rows = F, 
        use_raster = T, raster_quality = 10, 
        row_split = factor(c(rep('CLL common auto', length(CLL.common.auto)),
                      rep('CLL common mito', length(CLL.common.mito)),
                      rep('CLL pre-FCR auto', length(CLL.clone1.auto)),
                      rep('CLL pre-FCR mito', length(CLL.clone1.mito)),
                      rep('CLL post-RIC auto', length(CLL.clone2.auto)),
                      rep('CLL post-RIC mito', length(CLL.clone2.mito)),
                      rep('Donor mito', length(donor.mito))),
                      levels = c('CLL common auto', 'CLL common mito', 
                                 'CLL pre-FCR auto', 'CLL pre-FCR mito', 
                                 'CLL post-RIC auto', 'CLL post-RIC mito')),
        column_split = factor(c(rep('CLL pre-FCR', length(cells.1.CLL)),
                                rep('Immune cells pre-FCR', length(cells.1.immune)),
                                rep('CLL post-RIC', length(cells.2.recipient)), 
                                rep('Donor', length(cells.2.donor))), 
                              levels = c('CLL pre-FCR', 'CLL post-RIC', 'Immune cells pre-FCR', 'Donor')),
        column_title_gp = gpar(fontsize=8), 
        row_title_gp = gpar(fontsize=8),
        col = col_fun, border = T, 
        top_annotation = ha)
dev.off()

### plot selected mtDNA mutations on UMAP

df = as.data.frame(Embeddings(so.12, reduction = 'umap')[c(cells.1.CLL, cells.1.immune, cells.2.recipient, cells.2.donor),])
df = cbind(df, vafs.mito[rownames(df), ])

df$celltype = so.12$manual.cluster[rownames(df)]
df$sample = so.12$orig.ident[rownames(df)]

df = cbind(df, t(GetAssayData(so.12, slot = 'scale.data')[,rownames(df)]))

ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=celltype)) + 
  geom_point(size=0.5) +
  scale_color_manual(values = c('CLL' = 'orange', 'CD4 T cell' = 'lightblue', 'CD8 T cell' = 'darkblue', 'NK' = 'purple','Myeloid' = 'darkgreen',
                                'Cluster 6' = 'orange')) + 
  theme_classic() + 
  NoAxes() +
  theme(legend.position = 'none') 
ggsave('./figure_coevolution/CLL/figures/UMAPS/20230426_CLL4_UMAP_celltypes.png', width = 4, height = 4, dpi = 600)

ggplot(df[order(df$`3538G>A`),], aes(x=UMAP_1, y=UMAP_2, color=`3538G>A`)) + 
  geom_point(size=0.5) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'solar_rojos')) +
  theme_classic() + 
  NoAxes() +
  theme(legend.position = 'none') 
ggsave('./figure_coevolution/CLL/figures/UMAPS/20230426_CLL4_UMAP_3538G>A.png', width = 4, height = 4, dpi = 600)

ggplot(df[order(df$`16247A>G`),], aes(x=UMAP_1, y=UMAP_2, color=`16247A>G`)) + 
  geom_point(size=0.5) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'solar_rojos')) +
  theme_classic() + 
  NoAxes() +
  theme(legend.position = 'none') 
ggsave('./figure_coevolution/CLL/figures/UMAPS/20230426_CLL4_UMAP_16247A>G.png', width = 4, height = 4, dpi = 600)

ggplot(df[order(df$`16290C>T`),], aes(x=UMAP_1, y=UMAP_2, color=`16290C>T`)) + 
  geom_point(size=0.5) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'solar_rojos')) +
  theme_classic() + 
  NoAxes() +
  theme(legend.position = 'none') 
ggsave('./figure_coevolution/CLL/figures/UMAPS/20230426_CLL4_UMAP_16290C>T.png', width = 4, height = 4, dpi = 600)




ggplot(df[order(df$`1918G>A`),], aes(x=UMAP_1, y=UMAP_2, color=`1918G>A`)) + 
  geom_point(size=1) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'solar_rojos')) +
  theme_classic() + 
  NoAxes() +
  theme(legend.position = 'none') 
ggsave('./figure_coevolution/CLL/figures/UMAPS/20230426_CLL4_UMAP_1918G>A.png', width = 4, height = 4, dpi = 600)

ggplot(df[order(df$`5650G>A`),], aes(x=UMAP_1, y=UMAP_2, color=`5650G>A`)) + 
  geom_point(size=1) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'solar_rojos')) +
  theme_classic() + 
  NoAxes() +
  theme(legend.position = 'none') 
ggsave('./figure_coevolution/CLL/figures/UMAPS/20230426_CLL4_UMAP_5650G>A.png', width = 4, height = 4, dpi = 600)

ggplot(df[order(df$`786G>A`),], aes(x=UMAP_1, y=UMAP_2, color=`786G>A`)) + 
  geom_point(size=1) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'solar_rojos')) +
  theme_classic() + 
  NoAxes() +
  theme(legend.position = 'none') 
ggsave('./figure_coevolution/CLL/figures/UMAPS/20230426_CLL4_UMAP_786G>A.png', width = 4, height = 4, dpi = 600)

for (p in c('CD3', 'CD5', 'CD19', 'CD4', 'CD8', 'CD14', 'CD33', 'CD56')) {
  q=ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=df[,p])) + 
    geom_point(size=0.5) +
    scale_color_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius')) +
    theme_classic() + 
    NoAxes() +
    theme(legend.position = 'none') 
  ggsave(paste0('./figure_coevolution/CLL/figures/UMAPS/20230426_CLL4_UMAP_',p,'.png'), width = 4, height = 4, dpi = 600, plot = q)  
}