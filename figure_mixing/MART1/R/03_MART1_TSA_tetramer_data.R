library(ArchR)
library(ggplot2)
library(ggrepel)
library(parallel)
library(Seurat)

TSAB.so = readRDS(file='./data/mixing/MART1/objects/20220225_MART1_TSAB.so')

MART1.mito = loadArchRProject('./data/mixing/MART1/MART1.mito/')

# plot T cells
p=plotEmbedding(MART1.mito, name = 'Clusters', pal = c('C1' = 'grey', 'C2' = 'grey', 'C3' = 'black', 'C4' = 'grey', 'C5' = 'grey',
                                                     'C6' = 'black', 'C7' = 'black', 'C8' = 'black', 'C9' = 'black'))
plotPDF(list(p), ArchRProj = MART1.mito, name = 'Tcell.gate', addDOC = F)

# gate on T cells using chromatin accessibility profiles
Tcells = gsub(MART1.mito$cellNames[which(MART1.mito$Clusters %in% c('C3', 'C6', 'C7', 'C8', 'C9'))], 
                                  pattern = '-1', replacement = '')

# pull donor annotation from ArchR object
donor.df = data.frame(cell = MART1.mito$cellNames,
                      donor = MART1.mito$donor)
rownames(donor.df) = donor.df$cell

# include ADT data
df = as.data.frame(t(GetAssayData(TSAB.so, slot = 'scale.data')))
df$cluster = TSAB.so$seurat_clusters
df = df[Tcells,]
df$donor = donor.df[paste0(rownames(df),'-1'), 'donor']

# identify antigen-specific cells for MART1
cells = rownames(df)[which(df$SAV1 > 1.5 & df$SAV2 > 1.5 & df$SAV1 <5 & df$SAV2 < 5)]

ggplot() + geom_point(data=df[which(df$SAV1 < 5 & df$SAV2 < 5),], aes(x=SAV1, y=SAV2, color=Tcrb), size=0.5) +
  scale_color_gradientn(colours = c('white', BuenColors::jdb_palette(name = 'solar_rojos'))) + 
  scale_x_continuous('MART1 SAV1',limits = c(-2,5)) +
  scale_y_continuous('MART1 SAV2',limits = c(-2,5)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  theme_classic() +
  theme(axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/MART1/figures/20220226_MART1_SAV1_SAV2.svg', width = 3, height = 2.5)

# statistics
statistics.df = data.frame()
for (s in c('MART1_1', 'MART1_2', 'MART1_3')) {
  donor.1 =  rownames(df[which(df$donor == '108'),])
  donor.1 = donor.1[which(grepl(s, donor.1))]
  
  donor.2 =  rownames(df[which(df$donor == '5'),])
  donor.2 = donor.2[which(grepl(s, donor.2))]
  
  statistics.df = rbind(statistics.df, data.frame(sample = s, 
                                                  cells = length(which(grepl(s, MART1.mito$cellNames))),
                                                  Tcells = length(which(grepl(s, Tcells))),
                                                  donor.1 = length(donor.1),
                                                  donor.1.MART1 = length(which(cells %in% donor.1)),
                                                  donor.2 = length(donor.2),
                                                  donor.2.MART1 = length(which(cells %in% donor.2))))
}
statistics.df$donor.1.MART1.freq = statistics.df$donor.1.MART1 / statistics.df$cells
ggplot(statistics.df, aes(x=sample, y=100*donor.1.MART1.freq)) + 
  geom_line(aes(group=1)) + 
  geom_point(color='orange') + 
  scale_x_discrete('dilution',labels = c('1:3', '1:30', '1:300')) + 
  scale_y_log10('% MART1-specific T cells') +
  theme_classic() +
  theme(axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/MART1/figures/20220226_MART1_specific_cells_titration.svg', width = 2.5, height = 1.5)

boo = reshape2::melt(statistics.df[,c('sample', 'donor.1', 'donor.2', 'donor.1.MART1', 'donor.2.MART1')])
boo$sample.variable = paste0(boo$sample, '.', boo$variable)
boo$donor = c(rep('donor.1', 3), rep('donor.2', 3), rep('donor.1', 3), rep('donor.2', 3))
boo$type = c(rep('all', 6), rep('MART1', 6))
boo$donor.sample = factor(paste0(boo$donor, '.', boo$sample), 
                          levels = c('donor.1.MART1_1', 'donor.2.MART1_1', 'donor.1.MART1_2', 'donor.2.MART1_2',
                                     'donor.1.MART1_3', 'donor.2.MART1_3'))

ggplot(data=boo, aes(x=donor.sample, y=value, fill=variable, alpha=variable)) +
  scale_fill_manual(values = c('donor.1' = 'orange', 'donor.2' = 'blue', 'donor.1.MART1' = 'orange', 'donor.2.MART1' = 'blue')) +
  scale_alpha_manual(values = c('donor.1' = 0.3, 'donor.2' = 0.3, 'donor.1.MART1' = 1, 'donor.2.MART1' = 1)) + 
  geom_col() + 
  scale_y_sqrt('T cells', breaks = c(1,10,100,500, 1000,2000)) +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave('./figure_mixing/MART1/figures/20220227_MART1_cells.svg', width = 1.5, height = 1.5)

TSAB.so$manual.cluster = 'none'
TSAB.so$manual.cluster[which(TSAB.so$seurat_clusters %in% c('0', '1', '2', '4', '13', '19'))] = 'CD4+ T cells'
TSAB.so$manual.cluster[which(TSAB.so$seurat_clusters %in% c('3','8', '9', '10', '11'))] = 'CD8+ T cells'
TSAB.so$manual.cluster[which(TSAB.so$seurat_clusters %in% c('5', '6', '7', '12', '18'))] = 'Monocytes'
TSAB.so$manual.cluster[which(TSAB.so$seurat_clusters %in% c('14'))] = 'B cells'
TSAB.so$manual.cluster[which(TSAB.so$seurat_clusters %in% c('20'))] = 'CD8+ T cells MART1'
TSAB.so$manual.cluster[which(TSAB.so$seurat_clusters %in% c('22'))] = 'CD4+ T cells MART1'

Idents(TSAB.so) = 'manual.cluster'
boo = subset(TSAB.so, downsample=200)
boo$manual.cluster = factor(boo$manual.cluster, levels = c('B cells', 'Monocytes', 'CD4+ T cells', 'CD4+ T cells MART1', 'CD8+ T cells', 'CD8+ T cells MART1'))
DoHeatmap(boo, cells = colnames(boo)[which(boo$manual.cluster != 'none')], group.by = 'manual.cluster', 
          features = c('SAV1', 'SAV2', 'Tcrb','CD3', 'CD4', 'CD8', 'CD127', 
                       'CD45RO', 'CD45RA', 'CD62L', 'CCR7', 'CD14', 'CD33', 'CD19'), 
          group.colors = c('B cells' = 'red', 'Monocytes' = 'darkgreen', 
                           'CD4+ T cells' = RColorBrewer::brewer.pal(name = 'Blues', n=3)[2],
                           'CD4+ T cells MART1' = 'grey',
                           'CD8+ T cells' = RColorBrewer::brewer.pal(name = 'Blues', n=3)[3],
                           'CD8+ T cells MART1' = 'black'), size = 3) + 
  scale_fill_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_yes', type='continuous'), na.value = 'white') +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/MART1/figures/20220227_MART1_heatmap.svg', width = 5, height = 2.5)

cells = colnames(TSAB.so)[which(TSAB.so$manual.cluster %in% c('CD4+ T cells MART1', 'CD8+ T cells MART1'))]
cells = cells[paste0(cells, '-1') %in% MART1.mito$cellNames]
df$donor[which(rownames(df) %in% cells)]
