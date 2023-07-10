library(ArchR)
library(dplyr)
library(gplots)
library(parallel)
library(Seurat)
library(ComplexHeatmap)

setwd('/Users/liviuspenter/dfci/asap_seq/')
source('./analysis/Tcell.colors.R')
source('/Users/shaka87/dfci/scripts/20200328_mtscatac_seq.R')

AML.1026.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1026.mito/')
TSB.so = readRDS('./data/objects/20210601_TSB.so.rds')

AML.1026.cluster.colors = c('0' = BuenColors::jdb_palette(name = 'corona')[1],
                            '1' = BuenColors::jdb_palette(name = 'corona')[2],
                            '2' = BuenColors::jdb_palette(name = 'corona')[3],
                            '3' = BuenColors::jdb_palette(name = 'corona')[4],
                            '4' = BuenColors::jdb_palette(name = 'corona')[5],
                            '5' = BuenColors::jdb_palette(name = 'corona')[6],
                            '6' = BuenColors::jdb_palette(name = 'corona')[7])

### 1026 
s = 'AML1026'

# read mtDNA data
variants = read.table(paste0('./data/mtDNA/20210429_',s,'_high_confidence_variants.csv'), header = T)
combined.frequencies = readRDS(paste0('./data/mtDNA/20210429_',s,'_combined_mutation_frequencies.rds'))

# cells with mtDNA information 
cells = AML.1026.mito$cellNames[which(grepl(s,AML.1026.mito$Sample) & !AML.1026.mito$manual.clusters %in% c('T cell'))]
#cells = AML.1026.mito$cellNames[which(grepl(s,AML.1026.mito$Sample))]
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

chimerism.df = read.csv2('./data/mtDNA/20200302_AML1026_donor_recipient.csv', row.names = 1)

df = uwot::umap(t(combined.frequencies), n_threads = 12)
df = data.frame(df)
colnames(df) = c('UMAP1', 'UMAP2')
ggplot(df, aes(x=UMAP1, y=UMAP2)) + geom_point(size=0.5)

mtDNA.so = CreateSeuratObject(combined.frequencies)
mtDNA.so = FindVariableFeatures(mtDNA.so)
mtDNA.so = ScaleData(mtDNA.so, features = rownames(mtDNA.so))
mtDNA.so = RunPCA(mtDNA.so)
mtDNA.so = RunUMAP(mtDNA.so, dims = 1:10)
mtDNA.so = RunTSNE(mtDNA.so)
mtDNA.so = FindNeighbors(mtDNA.so)
mtDNA.so = FindClusters(mtDNA.so, resolution = 0.1)

# get rid of donor-derived cells
mtDNA.so = subset(mtDNA.so, seurat_clusters != '8')

mtDNA.markers <- FindAllMarkers(mtDNA.so, min.pct = 0.05, logfc.threshold = 0.05)
mtDNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(mtDNA.so, downsample = 100), features = top10$gene, slot = 'counts', disp.max = 0.1,
          group.colors = AML.1026.cluster.colors, label = F) + NoLegend() + 
  scale_fill_gradientn(colours = c('white', BuenColors::jdb_palette(name = 'solar_rojos')), na.value = 'black') +
  theme(axis.text = element_text('Arial', size=8, color='black'))
ggsave('./figures/heatmaps/20220304_AML1026_mtDNA_clusters.svg', width = 3, height = 5)
ggsave('./figures/heatmaps/20220304_AML1026_mtDNA_clusters.png', width = 3, height = 2.5, dpi=600)


# mtDNA clusters versus celltypes
AML.1026.mito$mito.cluster = mtDNA.so$seurat_clusters[AML.1026.mito$cellNames]
df = data.frame(bc = AML.1026.mito$cellNames,
                sample = AML.1026.mito$Sample.hash,
                manual.cluster = AML.1026.mito$manual.clusters,
                mito.cluster = AML.1026.mito$mito.cluster)
df$UMAP1 = getEmbedding(AML.1026.mito)[df$bc,1]
df$UMAP2 = getEmbedding(AML.1026.mito)[df$bc,2]

df2 <- df %>% filter(!is.na(mito.cluster)) %>% group_by(sample, manual.cluster, mito.cluster) %>% tally()
df2 <- df2 %>% group_by(sample, manual.cluster) %>% mutate(freq = n / sum(n))
df2$manual.cluster <- factor(df2$manual.cluster, levels = c('HSCT', 'GMP', 'Mono', 'erythroid', 'T cell'))
df2$mito.cluster <- factor(df2$mito.cluster, levels = names(AML.1026.cluster.colors))

for (s in c('AML1026_1', 'AML1026_2', 'AML1026_3', 'AML1026_4')) {
  p=ggplot(df2[which(df2$sample == s & df2$manual.cluster %in% c('HSCT', 'GMP', 'Mono', 'erythroid')),], 
           aes(x=manual.cluster, y=100*freq, fill=mito.cluster)) + 
    geom_col() +
    scale_fill_manual(values = AML.1026.cluster.colors) + 
    scale_x_discrete(labels = c('HSC', 'GMP', 'Mono', 'erythroid')) + 
    scale_y_continuous('% cells') + 
    theme_classic() + 
    theme(legend.position = 'none',
          axis.title = element_text('Arial', size=10, color='black'),
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  ggsave(paste0('./figures/plots/AML1026/20230331_mito_cluster_dynamics_', s, '.svg'), width=1.3, height = 2, plot = p)
}

for (s in c('AML1026_1', 'AML1026_2', 'AML1026_3', 'AML1026_4')) {
  p=ggplot(df[which(df$sample == s),], aes(x=UMAP1, y=UMAP2, color=mito.cluster)) + geom_point(size=0.5) + 
    scale_color_manual(values = AML.1026.cluster.colors) + 
    theme_classic() +
    NoAxes() + NoLegend()
  ggsave(paste0('./figures/umaps/AML1026/20230331_',s,'_UMAP_mitoclusters.png'), width = 3, height = 3, dpi = 600, plot = p)
}


# small plot for clinical annotation 
boo = data.frame(label = c('Screening', 'Decitabine', 'C1', 'EOT'),
                 blasts = c(18,17,3,18))
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
ggsave('./figures/plots/AML1026/20230331_clinical_course.svg', width = 1.3, height = 1.5)

boo = data.frame(bc = AML.1026.mito$cellNames, cluster = AML.1026.mito$manual.clusters, sample.hash = AML.1026.mito$Sample.hash)
rownames(boo) = boo$bc
mtDNA.so$manual.cluster = boo[colnames(mtDNA.so), 'cluster']
mtDNA.so$sample.hash = boo[colnames(mtDNA.so), 'sample.hash']
TSB.so = RenameCells(TSB.so, new.names = paste0(colnames(TSB.so), '-1'))
boo = TSB.so[,which(colnames(TSB.so) %in% colnames(mtDNA.so))]
mtDNA.so[['ADT']] = boo@assays$ADT

boo = as.data.frame(prop.table(table(Idents(mtDNA.so), mtDNA.so$sample.hash), margin = 2))

ggplot(boo[which(boo$Var2 != 'none'),], 
       aes(x=Var2, y=100*Freq, color=Var1)) + geom_line(aes(group=Var1)) +
  scale_x_discrete(labels = c('Screening', 'EOLN', 'C1', 'EOT')) + 
  scale_y_log10('% cells') + 
  scale_color_manual(values = AML.1026.cluster.colors) + 
  theme_classic() +
  theme(legend.position = 'none', 
        axis.title.x = element_blank(),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text('Arial', size=10, color='black', angle=90, hjust=1, vjust=0.5),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figures/plots/20220303_AML1026_mtDNA_clusters.svg', width = 1.2, height = 2)

DimPlot(mtDNA.so) +
  scale_color_manual(values = AML.1026.cluster.colors) +
  NoAxes() + NoLegend()
ggsave('./figures/umaps/20220304_AML1026_mtDNA_clusters.png', width = 3, height = 3, dpi = 600)

saveRDS(mtDNA.so, file = './data/objects/20200303_AML1026_mtDNA.so.rds')
