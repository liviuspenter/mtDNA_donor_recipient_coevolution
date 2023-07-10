# explore donor chimerism across time 

library(ArchR)
library(ggplot2)
library(ggrepel)
library(gplots)
library(ggsignif)
library(parallel)
library(Signac)
library(Seurat)
library(BuenColors)
library(dplyr)

source('./figure_IST/R/IST.sample.colors.R')

IST.asap.mito = loadArchRProject('./data/IST/IST.asap.mito/')
TSB.so = readRDS('./data/IST/objects/20220117_IST_TSB.so.rds')
TSB.so=RenameCells(TSB.so, new.names = paste0(colnames(TSB.so), '-1'))
TSB.so.mito = subset(TSB.so, cells=colnames(TSB.so)[which(colnames(TSB.so) %in% IST.asap.mito$cellNames)])

# add TSB information
TSB.so.mito <- NormalizeData(TSB.so.mito, assay = "ADT", normalization.method = 'CLR', margin = 2)
TSB.so.mito = ScaleData(TSB.so.mito)
TSB.so.mito = FindVariableFeatures(TSB.so.mito)
TSB.so.mito = RunPCA(TSB.so.mito)
TSB.so.mito = RunUMAP(TSB.so.mito, dims = 1:20)
TSB.so.mito = FindNeighbors(TSB.so.mito)
TSB.so.mito = FindClusters(TSB.so.mito)
saveRDS(file='./data/objects/20220118_IST_mito_TSB.so.rds', TSB.so.mito)

# annotate clusters in ArchR object
IST.asap.mito$manual.cluster = 'none'
IST.asap.mito$manual.cluster[which(IST.asap.mito$Clusters %in% c('C8','C9'))] = 'B cell'
IST.asap.mito$manual.cluster[which(IST.asap.mito$Clusters %in% c('C2'))] = 'HSC'
IST.asap.mito$manual.cluster[which(IST.asap.mito$Clusters %in% c('C1'))] = 'Erythroid'
IST.asap.mito$manual.cluster[which(IST.asap.mito$Clusters %in% c('C3','C4','C5', 'C6', 'C7'))] = 'Myeloid'
IST.asap.mito$manual.cluster[which(IST.asap.mito$Clusters %in% c('C12', 'C13'))] = 'NK'
IST.asap.mito$manual.cluster[which(IST.asap.mito$Clusters %in% c('C10', 'C11'))] = 'T'
saveArchRProject(IST.asap.mito)

# plot TSB expression
p=plotEmbedding(IST.asap.mito, name = 'CD19')
plotPDF(list(p), name = 'CD19', addDOC = F, ArchRProj = IST.asap.mito)
p=plotEmbedding(IST.asap.mito, name = 'CD138')
plotPDF(list(p), name = 'CD138', addDOC = F, ArchRProj = IST.asap.mito)
p=plotEmbedding(IST.asap.mito, name = 'CD14')
plotPDF(list(p), name = 'CD14', addDOC = F, ArchRProj = IST.asap.mito)
p=plotEmbedding(IST.asap.mito, name = 'CD16')
plotPDF(list(p), name = 'CD16', addDOC = F, ArchRProj = IST.asap.mito)
p=plotEmbedding(IST.asap.mito, name = 'CD33')
plotPDF(list(p), name = 'CD33', addDOC = F, ArchRProj = IST.asap.mito)
p=plotEmbedding(IST.asap.mito, name = 'CD3')
plotPDF(list(p), name = 'CD3', addDOC = F, ArchRProj = IST.asap.mito)
p=plotEmbedding(IST.asap.mito, name = 'CD4')
plotPDF(list(p), name = 'CD4', addDOC = F, ArchRProj = IST.asap.mito)
p=plotEmbedding(IST.asap.mito, name = 'CD8')
plotPDF(list(p), name = 'CD8', addDOC = F, ArchRProj = IST.asap.mito)
p=plotEmbedding(IST.asap.mito, name = 'CD56')
plotPDF(list(p), name = 'CD56', addDOC = F, ArchRProj = IST.asap.mito)
p=plotEmbedding(IST.asap.mito, name = 'CD117')
plotPDF(list(p), name = 'CD117', addDOC = F, ArchRProj = IST.asap.mito)

# plot donor-recipient
p=plotEmbedding(IST.asap.mito, name = 'individual', pal = c('donor' = 'orange', 'recipient' = 'purple', 'none' = 'grey'))
plotPDF(list(p), name = 'donor', addDOC = F, ArchRProj = IST.asap.mito)

p=plotEmbedding(IST.asap.mito, name = 'manual.cluster', pal = IST.cell.colors)
plotPDF(list(p), name = 'celltypes', addDOC = F, ArchRProj = IST.asap.mito)

# plot cell type kinetics
df = data.frame(sample = IST.asap.mito$Sample, 
                celltype = IST.asap.mito$manual.cluster)

df = df %>% group_by(sample, celltype) %>% summarize(n = n()) %>%
  mutate(freq = n / sum(n))
df$celltype = factor(df$celltype, levels = c('Myeloid', 'HSC', 'Erythroid', 'B cell','NK', 'T'))
ggplot(df, aes(x=sample, y=100*freq, fill=celltype)) + geom_col() + 
  scale_x_discrete(limits = c('IST1_1', 'IST1_2', 'IST2_1', 'IST2_2', 'IST3_1', 'IST3_2', 'IST4_1', 'IST4_2', 'IST5_1', 'IST5_2'),
                   labels = c('PB pre', 'PB post', 'BM pre', 'BM post', 'BM pre', 'BM post', 'PB pre', 'PB post', 'PB pre', 'PB post')) +
  scale_y_continuous('% cells') + 
  scale_fill_manual(values = IST.cell.colors) +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 
ggsave('./figure_IST/figures/plots/20220327_celltype_kinetics.svg', width = 2, height = 2)

df$timepoint = factor(ifelse(grepl('_2', df$sample), 'post', 'pre'), levels = c('pre', 'post')) 
df$group = stringr::str_split_fixed(df$sample, pattern = '_', n=2)[,1]
ggplot(df[which(df$celltype == 'Myeloid' & df$group != 'IST2'),], aes(x=timepoint, y=100*freq)) + 
  geom_point(size=0.5, color=IST.cell.colors['Myeloid']) + 
  geom_signif(comparisons = list(c('pre', 'post')), textsize = 2, test = 't.test', test.args = c('alternative' = 'less')) + 
  geom_line(aes(group=group), color=IST.cell.colors['Myeloid']) +
  scale_x_discrete(labels = c('Pre', 'Post')) + 
  scale_y_continuous('% cells') + 
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 
ggsave('./figure_IST/figures/plots/20220327_celltype_kinetics_myeloid.svg', width = 1, height = 1.5)

ggplot(df[which(df$celltype == 'T'& df$group != 'IST2'),], aes(x=timepoint, y=100*freq)) + 
  geom_point(size=0.5, color=IST.cell.colors['T']) + 
  geom_signif(comparisons = list(c('pre', 'post')), textsize = 2, test = 't.test', test.args = c('alternative' = 'greater')) + 
  geom_line(aes(group=group), color=IST.cell.colors['T']) +
  scale_x_discrete(labels = c('Pre', 'Post')) + 
  scale_y_continuous('% cells') + 
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 
ggsave('./figure_IST/figures/plots/20220327_celltype_kinetics_T.svg', width = 1, height = 1.5)

# plot cell type-specific donor chimerism
chimerism.df = data.frame()
for (mc in unique(IST.asap.mito$manual.cluster)) {
  for (sample in unique(IST.asap.mito$Sample)) {
    donor.cells = length(which(IST.asap.mito$manual.cluster == mc & IST.asap.mito$Sample == sample & IST.asap.mito$individual == 'donor'))
    cells = length(which(IST.asap.mito$manual.cluster == mc & IST.asap.mito$Sample == sample & IST.asap.mito$individual %in% c('donor', 'recipient')))
    chimerism.df = rbind(chimerism.df, data.frame(cluster = mc, Sample = sample, donor.cells = donor.cells, cells = cells))
  }
}
chimerism.df$chimerism.freq = chimerism.df$donor.cells / chimerism.df$cells
chimerism.df$timepoint = stringr::str_split_fixed(chimerism.df$Sample, pattern = '_', n=2)[,2]
chimerism.df$patient = stringr::str_split_fixed(chimerism.df$Sample, pattern = '_', n=2)[,1]
chimerism.df$cluster = factor(chimerism.df$cluster, levels = c('HSC', 'Myeloid', 'T', 'NK', 'Erythroid', 'B cell'))
ggplot(chimerism.df[which(chimerism.df$cells > 10),], aes(x=timepoint, y=100*chimerism.freq, color=patient)) + geom_point(size=0.5) +
  geom_line(aes(group=patient)) + 
  scale_x_discrete(breaks = c(1,2), labels = c('Pre', 'Post')) + 
  scale_y_continuous('% donor chimerism') + 
  scale_color_manual(values = IST.sample.colors) + 
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 10)) +
  facet_wrap(~cluster)
ggsave('./figure_IST/figures/plots/20220118_IST_chimerism_kinetics.svg', width = 2.5, height = 2.5)

chimerism.df[which(chimerism.df$cells < 5), 'chimerism.freq'] = NA

for (celltype in unique(chimerism.df$cluster)) {
  p=ggplot(chimerism.df[which(chimerism.df$cluster == celltype),], aes(x=Sample, y=100*chimerism.freq)) + 
    geom_col(aes(alpha=as.character(timepoint)), fill = IST.cell.colors[celltype]) +
    scale_alpha_manual(values = c('1' = 0.5, '2' = 1)) + 
    scale_x_discrete(labels = c('PB pre', 'PB post', 'BM pre', 'BM post', 'BM pre', 'BM post', 
                                'PB pre', 'PB post', 'PB pre', 'PB post')) +
    scale_y_continuous('% donor') + 
    theme_classic() +
    theme(legend.position = 'none', 
          axis.text = element_text('Arial', size=8, color='black'),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          axis.title = element_text('Arial', size=8, color='black'),
          axis.title.x = element_blank()) 
  ggsave(paste0('./figure_IST/figures/plots/20230424_IST_donor_chimerism_', celltype, '.svg'), width = 1.5, height = 1.5, plot = p)
}

# hone in on TNK cells

IST.asap.mito.TNK = subsetArchRProject(ArchRProj = IST.asap.mito, 
                                   cells = IST.asap.mito$cellNames[which(IST.asap.mito$manual.cluster %in% c('T','NK'))],
                                   outputDirectory = './data/IST/IST.asap.mito.TNK/', threads = 12, force = T)
IST.asap.mito.TNK = addClusters(input = IST.asap.mito.TNK, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 0.3, force = T)
IST.asap.mito.TNK = addUMAP(ArchRProj = IST.asap.mito.TNK, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
IST.asap.mito.TNK = addImputeWeights(IST.asap.mito.TNK)

IST.asap.mito.TNK = addHarmony(IST.asap.mito.TNK, reducedDims = 'IterativeLSI', name = 'Harmony', groupBy = 'Sample')
IST.asap.mito.TNK = addUMAP(IST.asap.mito.TNK, reducedDims = 'Harmony', name = 'UMAP.harmony', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
IST.asap.mito.TNK = addClusters(input = IST.asap.mito.TNK, reducedDims = 'Harmony', method = 'Seurat', name = 'Clusters.harmony', resolution = 1, force = T)
IST.asap.mito.TNK$manual.cluster = 'NK'
IST.asap.mito.TNK$manual.cluster[which(IST.asap.mito.TNK$Clusters.harmony %in% c('C9','C10', 'C11', 'C12'))] = 'CD4'
IST.asap.mito.TNK$manual.cluster[which(IST.asap.mito.TNK$Clusters.harmony %in% c('C2'))] = 'CD8'
saveArchRProject(IST.asap.mito.TNK)
IST.asap.mito.TNK = loadArchRProject('./data/IST/IST.asap.mito.TNK/')

TSB.so.mito.TNK = subset(TSB.so, cells=colnames(TSB.so)[which(colnames(TSB.so) %in% IST.asap.mito.TNK$cellNames)])
TSB.so.mito.TNK <- NormalizeData(TSB.so.mito.TNK, assay = "ADT", normalization.method = 'CLR', margin = 2)
TSB.so.mito.TNK = ScaleData(TSB.so.mito.TNK)
TSB.so.mito.TNK = FindVariableFeatures(TSB.so.mito.TNK)
TSB.so.mito.TNK = RunPCA(TSB.so.mito.TNK)
TSB.so.mito.TNK = RunUMAP(TSB.so.mito.TNK, dims = 1:20)
TSB.so.mito.TNK = FindNeighbors(TSB.so.mito.TNK)
TSB.so.mito.TNK = FindClusters(TSB.so.mito.TNK)
saveRDS(file='./data/objects/20220118_IST_mito_TNK_TSB.so.rds', TSB.so.mito.TNK)

# plot donor chimerism only within TNK cells
chimerism.df.TNK = data.frame()
for (mc in unique(IST.asap.mito.TNK$manual.cluster)) {
  for (sample in unique(IST.asap.mito.TNK$Sample)) {
    donor.cells = length(which(IST.asap.mito.TNK$manual.cluster == mc & IST.asap.mito.TNK$Sample == sample & IST.asap.mito.TNK$individual == 'donor'))
    cells = length(which(IST.asap.mito.TNK$manual.cluster == mc & IST.asap.mito.TNK$Sample == sample & IST.asap.mito.TNK$individual %in% c('donor', 'recipient')))
    chimerism.df.TNK = rbind(chimerism.df.TNK, data.frame(cluster = mc, Sample = sample, donor.cells = donor.cells, cells = cells))
  }
}
chimerism.df.TNK$chimerism.freq = chimerism.df.TNK$donor.cells / chimerism.df.TNK$cells
chimerism.df.TNK$chimerism.freq[which(chimerism.df.TNK$cells < 5)] = NA
chimerism.df.TNK$timepoint = stringr::str_split_fixed(chimerism.df.TNK$Sample, pattern = '_', n=2)[,2]
chimerism.df.TNK$patient = stringr::str_split_fixed(chimerism.df.TNK$Sample, pattern = '_', n=2)[,1]
chimerism.df.TNK$cluster = factor(chimerism.df.TNK$cluster, levels = c('CD4', 'CD8', 'NK'))
ggplot(chimerism.df.TNK, aes(x=timepoint, y=100*chimerism.freq, color=patient)) + geom_point(size=0.5) +
  geom_line(aes(group=patient)) + 
  scale_x_discrete(breaks = c(1,2), labels = c('Pre', 'Post')) + 
  scale_y_continuous('% donor chimerism') + 
  scale_color_manual(values = IST.sample.colors) + 
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 10)) +
  facet_wrap(~cluster)
ggsave('./figure_IST/figures/plots/20220118_IST_chimerism_kinetics_TNK.svg', width = 2, height = 2.5)

TNK.dynamics = data.frame()
for (sample in unique(IST.asap.mito.TNK$Sample)) {
  cells.sample = length(which(IST.asap.mito.TNK$Sample == sample))
  for (mc in unique(IST.asap.mito.TNK$manual.cluster)) {
    cells.cluster = length(which(IST.asap.mito.TNK$manual.cluster == mc & IST.asap.mito.TNK$Sample == sample))
    TNK.dynamics = rbind(TNK.dynamics, data.frame(sample = sample, manual.cluster = mc, cells.cluster = cells.cluster, 
                                                  cells.cluster.freq = cells.cluster / cells.sample))
  }
}
TNK.dynamics$timepoint = stringr::str_split_fixed(TNK.dynamics$sample, pattern = '_', n=2)[,2]
TNK.dynamics$patient = stringr::str_split_fixed(TNK.dynamics$sample, pattern = '_', n=2)[,1]

ggplot(TNK.dynamics, aes(x=timepoint, y=100*cells.cluster.freq, color=patient)) + geom_point(size=0.5) +
  geom_line(aes(group=patient)) + 
  scale_x_discrete(breaks = c(1,2), labels = c('Pre', 'Post')) + 
  scale_y_continuous('% T/NK cells') + 
  scale_color_manual(values = IST.sample.colors) + 
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 10)) +
  facet_wrap(~manual.cluster)
ggsave('./figure_IST/figures/plots/20220118_IST_kinetics_TNK.svg', width = 2, height = 2.5)

