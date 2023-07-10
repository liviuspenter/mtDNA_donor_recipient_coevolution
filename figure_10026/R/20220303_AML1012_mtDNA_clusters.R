library(ArchR)
library(dplyr)
library(gplots)
library(parallel)
library(Seurat)
library(ComplexHeatmap)

setwd('/Users/liviuspenter/dfci/asap_seq/')
source('./analysis/Tcell.colors.R')
source('/Users/shaka87/dfci/scripts/20200328_mtscatac_seq.R')

AML.1012.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1012.mito/')
TSB.so = readRDS('./data/objects/20210601_TSB.so.rds')

AML.1012.cluster.colors = c('0' = BuenColors::jdb_palette(name = 'corona')[1],
                            '1' = BuenColors::jdb_palette(name = 'corona')[2],
                            '2' = BuenColors::jdb_palette(name = 'corona')[3],
                            '3' = BuenColors::jdb_palette(name = 'corona')[4])

### 1012
s = 'AML1012'

# read mtDNA data
variants = read.table(paste0('./data/mtDNA/20210429_',s,'_high_confidence_variants.csv'), header = T)
combined.frequencies = readRDS(paste0('./data/mtDNA/20210429_',s,'_combined_mutation_frequencies.rds'))

# cells with mtDNA information 
cells = AML.1012.mito$cellNames[which(grepl(s,AML.1012.mito$Sample) & !AML.1012.mito$manual.clusters %in% c('CD4', 'CD8') )]
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

mtDNA.so = CreateSeuratObject(combined.frequencies)
mtDNA.so = FindVariableFeatures(mtDNA.so)
mtDNA.so = ScaleData(mtDNA.so, features = rownames(mtDNA.so))
mtDNA.so = RunPCA(mtDNA.so)
mtDNA.so = RunUMAP(mtDNA.so, dims = 1:20)
mtDNA.so = RunTSNE(mtDNA.so)
mtDNA.so = FindNeighbors(mtDNA.so)
mtDNA.so = FindClusters(mtDNA.so, resolution = 0.1)

boo = data.frame(bc = AML.1012.mito$cellNames, cluster = AML.1012.mito$manual.clusters, sample.hash = AML.1012.mito$Sample.hash)
rownames(boo) = boo$bc
mtDNA.so$manual.cluster = boo[colnames(mtDNA.so), 'cluster']
mtDNA.so$sample.hash = boo[colnames(mtDNA.so), 'sample.hash']
TSB.so = RenameCells(TSB.so, new.names = paste0(colnames(TSB.so), '-1'))
boo = TSB.so[,which(colnames(TSB.so) %in% colnames(mtDNA.so))]
mtDNA.so[['ADT']] = boo@assays$ADT

mtDNA.markers <- FindAllMarkers(mtDNA.so, min.pct = 0.05, logfc.threshold = 0.05)
mtDNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(mtDNA.so, downsample = 100), features = top10$gene, slot = 'counts', disp.max = 0.1, group.colors =  
            AML.1012.cluster.colors,
          label = F) + NoLegend() + 
  scale_fill_gradientn(colours = c('white', BuenColors::jdb_palette(name = 'solar_rojos')), na.value = 'black') +
  theme(axis.text = element_text('Arial', size=8, color='black'))
ggsave('./figures/heatmaps/20220304_AML1012_mtDNA_clusters.svg', width = 3, height = 2.5)
ggsave('./figures/heatmaps/20220304_AML1012_mtDNA_clusters.png', width = 2, height = 2)

# mtDNA clusters versus celltypes
AML.1012.mito$mito.cluster = mtDNA.so$seurat_clusters[AML.1012.mito$cellNames]
df = data.frame(bc = AML.1012.mito$cellNames,
                sample = AML.1012.mito$Sample.hash,
                manual.cluster = AML.1012.mito$manual.clusters,
                mito.cluster = AML.1012.mito$mito.cluster)
df$UMAP1 = getEmbedding(AML.1012.mito)[df$bc,1]
df$UMAP2 = getEmbedding(AML.1012.mito)[df$bc,2]

df2 <- df %>% group_by(sample, manual.cluster, mito.cluster) %>% tally()
df2 <- df2 %>% group_by(sample, manual.cluster) %>% mutate(freq = n / sum(n))
df2$manual.cluster <- factor(df2$manual.cluster, levels = c('HSCT', 'GMP', 'Mono', 'erythroid', 'CD4', 'CD8'))
df2$mito.cluster <- factor(df2$mito.cluster, levels = names(AML.1012.cluster.colors))

for (s in c('AML1012_1', 'AML1012_2', 'AML1012_3', 'AML1012_4')) {
  p=ggplot(df2[which(df2$sample == s & df2$manual.cluster %in% c('HSCT', 'GMP', 'Mono', 'erythroid')),], 
           aes(x=manual.cluster, y=100*freq, fill=mito.cluster)) + 
    geom_col() +
    scale_fill_manual(values = AML.1012.cluster.colors) + 
    scale_x_discrete(labels = c('HSC', 'GMP', 'Mono', 'erythroid')) + 
    scale_y_continuous('% cells') + 
    theme_classic() + 
    theme(legend.position = 'none',
          axis.title = element_text('Arial', size=10, color='black'),
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  ggsave(paste0('./figures/plots/AML1012/20230331_mito_cluster_dynamics_', s, '.svg'), width=1.3, height = 2, plot = p)
}

for (s in c('AML1012_1', 'AML1012_2', 'AML1012_3', 'AML1012_4')) {
  p=ggplot(df[which(df$sample == s),], aes(x=UMAP1, y=UMAP2, color=mito.cluster)) + geom_point(size=0.5) + 
    scale_color_manual(values = AML.1012.cluster.colors) + 
    theme_classic() +
    NoAxes() + NoLegend()
  ggsave(paste0('./figures/umaps/AML1012/20230331_',s,'_UMAP_mitoclusters.png'), width = 3, height = 3, dpi = 600, plot = p)
}

boo = as.data.frame(prop.table(table(Idents(mtDNA.so), mtDNA.so$sample.hash), margin = 2))

ggplot(boo[which(!boo$Var2 %in% c('none', 'AML1012_3')),], 
       aes(x=Var2, y=100*Freq, color=Var1)) + geom_line(aes(group=Var1)) +
  scale_x_discrete(labels = c('Screening', 'EOLN', 'EOT')) + 
  scale_y_log10('% cells') + 
  scale_color_manual(values = AML.1012.cluster.colors) + 
  theme_classic() +
  theme(legend.position = 'none', 
        axis.title.x = element_blank(),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text('Arial', size=10, color='black', angle=90, hjust=1, vjust=0.5),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figures/plots/20220303_AML1012_mtDNA_clusters.svg', width = 1.1, height = 2)

DimPlot(mtDNA.so) +
  scale_color_manual(values = AML.1012.cluster.colors) + 
  NoAxes() + NoLegend()
ggsave('./figures/umaps/20220304_AML1012_mtDNA_clusters.png', width = 3, height = 3, dpi = 600)

saveRDS(mtDNA.so, file = './data/objects/20200303_AML1012_mtDNA.so.rds')
saveArchRProject(AML.1012.mito)
