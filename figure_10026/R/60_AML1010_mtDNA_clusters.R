# discovery and tracking of AML subclones defined by mtDNA mutations in AML1010

library(ArchR)
library(dplyr)
library(gplots)
library(parallel)
library(Seurat)
library(ComplexHeatmap)

source('./figure_10026/R/Tcell.colors.R')
source('./R/20200328_mtscatac_seq.R')

AML.1010.cluster.colors = c('0' = BuenColors::jdb_palette(name = 'corona')[1],
                            '1' = BuenColors::jdb_palette(name = 'corona')[2],
                            '2' = BuenColors::jdb_palette(name = 'corona')[3],
                            '3' = BuenColors::jdb_palette(name = 'corona')[4],
                            '4' = BuenColors::jdb_palette(name = 'corona')[5],
                            '5' = BuenColors::jdb_palette(name = 'corona')[6],
                            '6' = BuenColors::jdb_palette(name = 'corona')[7],
                            '7' = BuenColors::jdb_palette(name = 'corona')[8],
                            '8' = BuenColors::jdb_palette(name = 'corona')[9],
                            '9' = BuenColors::jdb_palette(name = 'corona')[10],
                            '10' = BuenColors::jdb_palette(name = 'corona')[11],
                            '11' = BuenColors::jdb_palette(name = 'corona')[12],
                            '12' = BuenColors::jdb_palette(name = 'corona')[13])
AML.1010.cluster.colors = factor(AML.1010.cluster.colors, levels = AML.1010.cluster.colors)

AML.1010.mito = loadArchRProject('./data/10026/AML.1010.mito/')
TSB.so = readRDS('./data/10026/objects/20210601_TSB.so.rds')
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
cells = AML.1010.mito$cellNames[which(grepl(s,AML.1010.mito$Sample) & !AML.1010.mito$manual.clusters %in% c('CD4', 'CD8'))]
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

chimerism.df = read.csv2('./data/10026/mtDNA/20200302_AML1010_donor_recipient.csv', row.names = 1)

recipient.non.T.cells = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters %in% c('HSCT', 'GMP', 'Mono', 'erythroid'))]

mtDNA.so = CreateSeuratObject(combined.frequencies)
mtDNA.so = FindVariableFeatures(mtDNA.so)
mtDNA.so = ScaleData(mtDNA.so, features = rownames(mtDNA.so))
mtDNA.so = RunPCA(mtDNA.so)
mtDNA.so = RunUMAP(mtDNA.so, dims = 1:10)
mtDNA.so = RunTSNE(mtDNA.so)
mtDNA.so = FindNeighbors(mtDNA.so)
mtDNA.so = FindClusters(mtDNA.so, resolution = 0.1)

mtDNA.markers <- FindAllMarkers(mtDNA.so, min.pct = 0.05, logfc.threshold = 0.05)
mtDNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(mtDNA.so, downsample = 100), features = top10$gene, slot = 'counts', disp.max = 0.1, group.colors =  
            AML.1010.cluster.colors, label = F, raster = T) + #NoLegend() + 
  scale_fill_gradientn(colours = c('white', BuenColors::jdb_palette(name = 'solar_rojos')), na.value = 'black') +
  theme(axis.text = element_text('Arial', size=8, color='black'))
ggsave('./figure_10026/figures/heatmaps/20220304_AML1010_mtDNA_clusters.svg', width = 3, height = 5)
ggsave('./figure_10026/figures/heatmaps/20220304_AML1010_mtDNA_clusters.png', width = 3, height = 2.5, dpi=600)

# mtDNA clusters versus celltypes
AML.1010.mito$mito.cluster = mtDNA.so$seurat_clusters[AML.1010.mito$cellNames]
df = data.frame(bc = AML.1010.mito$cellNames,
                sample = AML.1010.mito$Sample.hash,
                manual.cluster = AML.1010.mito$manual.clusters,
                mito.cluster = AML.1010.mito$mito.cluster)
df$UMAP1 = getEmbedding(AML.1010.mito)[df$bc,1]
df$UMAP2 = getEmbedding(AML.1010.mito)[df$bc,2]

df2 <- df %>% group_by(sample, manual.cluster, mito.cluster) %>% tally()
df2 <- df2 %>% group_by(sample, manual.cluster) %>% mutate(freq = n / sum(n))
df2$manual.cluster <- factor(df2$manual.cluster, levels = c('HSCT', 'GMP', 'Mono', 'erythroid', 'CD4', 'CD8'))
df2$mito.cluster <- factor(df2$mito.cluster, levels = names(AML.1010.cluster.colors))

for (s in c('AML1010_1', 'AML1010_2', 'AML1010_3', 'AML1010_4', 'AML1010_5')) {
  p=ggplot(df2[which(df2$sample == s & df2$manual.cluster %in% c('HSCT', 'GMP', 'Mono', 'erythroid')),], 
         aes(x=manual.cluster, y=100*freq, fill=mito.cluster)) + 
    geom_col() +
    #scale_fill_manual(values = AML.1010.cluster.colors) + 
    scale_fill_manual(values = BuenColors::jdb_palette(name = 'corona',n=13)) + 
    scale_x_discrete(labels = c('HSC', 'GMP', 'Mono', 'erythroid')) + 
    scale_y_continuous('% cells') + 
    theme_classic() + 
    theme(legend.position = 'none',
          axis.title = element_text('Arial', size=10, color='black'),
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  ggsave(paste0('./figure_10026/figures/plots/AML1010/20230331_mito_cluster_dynamics_', s, '.svg'), width=1.3, height = 2, plot = p)
}

for (s in c('AML1010_1', 'AML1010_2', 'AML1010_3', 'AML1010_4', 'AML1010_5')) {
  p=ggplot(df[which(df$sample == s),], aes(x=UMAP1, y=UMAP2, color=mito.cluster)) + geom_point(size=0.5) + 
    scale_color_manual(values = BuenColors::jdb_palette(name = 'corona',n=13)) + 
    theme_classic() +
    NoAxes() + NoLegend()
  ggsave(paste0('./figure_10026/figures/umaps/AML1010/20230331_',s,'_UMAP_mitoclusters.png'), width = 3, height = 3, dpi = 600, plot = p)
}

# small plot for clinical annotation 
boo = data.frame(label = c('Screening', 'Decitabine', 'C1', 'C4', 'C10'),
                 blasts = c(37,21,34,38,73))
boo$label = factor(boo$label, levels = boo$label)
ggplot(boo, aes(x=label, y=blasts)) + geom_line(group=1) + 
  geom_point(size=0.5) +
  scale_y_continuous('% blasts',limits = c(0,100)) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./figure_10026/figures/plots/AML1010/20230331_clinical_course.svg', width = 1.3, height = 1.5)

boo = data.frame(bc = AML.1010.mito$cellNames, cluster = AML.1010.mito$manual.clusters, sample.hash = AML.1010.mito$Sample.hash)
rownames(boo) = boo$bc
mtDNA.so$manual.cluster = boo[colnames(mtDNA.so), 'cluster']
mtDNA.so$sample.hash = boo[colnames(mtDNA.so), 'sample.hash']
TSB.so = RenameCells(TSB.so, new.names = paste0(colnames(TSB.so), '-1'))
boo = TSB.so[,which(colnames(TSB.so) %in% colnames(mtDNA.so))]
mtDNA.so[['ADT']] = boo@assays$ADT

boo = as.data.frame(prop.table(table(Idents(mtDNA.so), mtDNA.so$sample.hash), margin = 2))

ggplot(boo[which(boo$Var2 != 'none'),], 
       aes(x=Var2, y=100*Freq, color=Var1)) + geom_line(aes(group=Var1)) +
  scale_x_discrete(labels = c('Screening', 'Decitabine', 'C1', 'C4', 'C10')) + 
  scale_y_log10('% cells') + 
  scale_color_manual(values = BuenColors::jdb_palette(name = 'corona',n=13)) + 
  theme_classic() +
  theme(legend.position = 'none', 
        axis.title.x = element_blank(),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text('Arial', size=10, color='black', angle=90, hjust=1, vjust=0.5),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_10026/figures/plots/20220303_AML1010_mtDNA_clusters.svg', width = 1.3, height = 2)

DimPlot(mtDNA.so) +
  scale_color_manual(values = BuenColors::jdb_palette(name = 'corona',n=13)) + 
  NoAxes() + NoLegend()
ggsave('./figure_10026/figures/umaps/20220304_AML1010_mtDNA_clusters.png', width = 3, height = 3, dpi = 600)

mtDNA.so$meta.cluster = 'Cluster_1'
mtDNA.so$meta.cluster[which(mtDNA.so$seurat_clusters %in% c('0','2'))] = 'Cluster_2'
Idents(mtDNA.so) = 'meta.cluster'
boo = prop.table(table(Idents(mtDNA.so), mtDNA.so$sample.hash), margin = 2)

ggplot(as.data.frame(t(boo[,c('AML1010_1', 'AML1010_2', 'AML1010_3', 'AML1010_4', 'AML1010_5')])), 
       aes(x=Var1, y=100*Freq, color=Var2)) + geom_line(aes(group=Var2)) +
  scale_x_discrete(labels = c('Screening', 'EOLN', 'C1', 'C4', 'C10')) + 
  scale_y_continuous('% cells') + 
  scale_color_manual(values = c('Cluster_1' = 'orange', 'Cluster_2' = 'red', 'T cell' = 'blue')) + 
  theme_classic() +
  theme(legend.position = 'none', 
        axis.title.x = element_blank(),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text('Arial', size=10, color='black', angle=90, hjust=1, vjust=0.5),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_10026/figures/plots/20220303_AML1010_mtDNA_clusters.svg', width = 1.3, height = 2)

df = data.frame()
for (s in c('AML1010_1', 'AML1010_2', 'AML1010_3', 'AML1010_4', 'AML1010_5')) {
  boo = subset(mtDNA.so, sample.hash == s & meta.cluster %in% c('Cluster_1', 'Cluster_2'))
  Idents(boo) = 'manual.cluster'
  boo2 = data.frame(prop.table(table(Idents(boo), boo$meta.cluster), margin = 1))
  boo2$sample = s
  df = rbind(df, boo2)
}
df$Freq = as.numeric(df$Freq)
boo = reshape2::dcast(df, sample ~ Var1 + Var2, value.var = 'Freq')
col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = BuenColors::jdb_palette(name = 'brewer_celsius'))
ha = HeatmapAnnotation(Cluster = c(rep('Cluster 1',4), rep('Cluster 2',4)), 
                       col = list(Cluster = c('Cluster 1' = 'orange', 'Cluster 2' = 'red')), simple_anno_size = unit(5, 'pt'), border=T)
svglite::svglite('./figure_10026/figures/heatmaps/20220303_AML1010_mtDNA_clusters.svg', width = 3, height = 2)
Heatmap(as.matrix(boo[,c('HSCT_Cluster_1',  'GMP_Cluster_1', 'Mono_Cluster_1', 
               'erythroid_Cluster_1', 'HSCT_Cluster_2', 'GMP_Cluster_2', 'Mono_Cluster_2', 'erythroid_Cluster_2')]),
        cluster_rows = F, cluster_columns = F, col = col_fun, column_split = c(rep('Cluster 1',4), rep('Cluster 2',4)), border = T, 
        column_labels = c('HSC', 'GMP', 'Mono', 'erythroid', 'HSC', 'GMP', 'Mono', 'erythroid'), top_annotation = ha, 
        column_names_gp = gpar(fontsize=10), column_title_gp = gpar(fontsize=10))
dev.off()

saveRDS(mtDNA.so, file = './data/10026/objects/20220303_AML1010_mtDNA.so.rds')
saveArchRProject(AML.1010.mito)