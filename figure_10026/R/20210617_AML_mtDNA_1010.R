library(ArchR)
library(gplots)
library(Seurat)
library(ComplexHeatmap)

setwd('/Users/liviuspenter/dfci/asap_seq/')
source('./analysis/Tcell.colors.R')

AML.1010.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1010.mito/')
TSB.so = readRDS('./data/objects/20210601_TSB.so.rds')

### 1010 
s = 'AML1010'
donor.variants = gtools::mixedsort(c('195T>C', '150C>T', '6146A>G', '1811A>G', '10907T>C', '6047A>G', '5999T>C', '14866C>T', '11009T>C', '11467A>G', '9070T>G', '15693T>C', '12308A>G',
                   '12372G>A', '4646T>C', '14620C>T', '15530T>C', '4811A>G', '499G>A', '16356T>C', '11332C>T', '16179C>T'))
recipient.variants = gtools::mixedsort(c('16126T>C', '4216T>C', '10084T>C', '489T>C', '462C>T', '11251A>G', '12612A>G', '15452C>A', '3010G>A', '14798T>C', '13708G>A', 
                       '16069C>T', '16319G>A', '295C>T', '55T>C', '56A>G', '185G>A', '228G>A'))

# read mtDNA data
variants = read.table(paste0('./data/mtDNA/20210429_',s,'_high_confidence_variants.csv'), header = T)
combined.frequencies = readRDS(paste0('./data/mtDNA/20210429_',s,'_combined_mutation_frequencies.rds'))
combined.frequencies = combined.frequencies[-which(rownames(combined.frequencies) == '310T>C'),]

# cells with mtDNA information
cells = AML.1010.mito$cellNames[which(grepl(s,AML.1010.mito$Sample))]
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

# filter uninformative mtDNA mutations
#combined.frequencies = combined.frequencies[-which(rowSums(combined.frequencies) == 0),]
combined.frequencies = combined.frequencies[-which(rowMeans(combined.frequencies) < 0.01),]
#combined.frequencies = combined.frequencies[variants$variant[which(variants$n_cells_conf_detected > 50)],]
#combined.frequencies = combined.frequencies[order(rowVars(as.matrix(combined.frequencies)), decreasing = T)[1:30],]

# order by cell phenotype
cells.ordered = c()
for (cluster in names(AML.1010.mito.colors)) {
  boo = AML.1010.mito$cellNames[which(AML.1010.mito$manual.clusters == cluster & AML.1010.mito$cellNames %in% cells)]
  names(boo) = rep(cluster, length(boo))
  cells.ordered = c(cells.ordered, boo)
}

chimerism.df = data.frame(donor.variants = colMeans(combined.frequencies[donor.variants,cells.ordered]),
                          recipient.variants = colMeans(combined.frequencies[recipient.variants,cells.ordered]))
chimerism.df$individual = 'none'
chimerism.df$individual[which(chimerism.df$donor.variants > 0.8 & chimerism.df$recipient.variants < 0.2)] = 'donor'
chimerism.df$individual[which(chimerism.df$donor.variants < 0.8 & chimerism.df$recipient.variants > 0.8)] = 'recipient'

p=ggplot(chimerism.df, aes(x=100*donor.variants, y=100*recipient.variants, color=individual)) + 
  ggrastr::rasterize(geom_point(size=0.5), dpi=600) +
  scale_x_continuous('% recipient variants') + 
  scale_y_continuous('% donor variants') + 
  scale_color_manual(values = c('none' = 'grey', 'donor' = 'purple', 'recipient' = 'orange')) + 
  geom_hline(yintercept = c(20,80)) + 
  geom_vline(xintercept = c(20,80)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave(paste0('./figures/deconvolution/20220302_AML1010_mtDNA_variants.svg'), width = 1.5, height = 1.5, plot = p)

df = getEmbedding(AML.1010.mito)
colnames(df) = c('UMAP1', 'UMAP2')
df$origin = chimerism.df[rownames(df),'individual']
df[AML.1010.mito$cellNames,'celltype'] = AML.1010.mito$manual.clusters
df[AML.1010.mito$cellNames,'Sample'] = AML.1010.mito$Sample.hash

ggplot() + 
  geom_point(data=df[which(df$origin == 'recipient'),], aes(x=UMAP1, y=UMAP2), color='lightgrey', size=0.5) +
  geom_point(data=df[which(df$origin == 'donor' & df$celltype %in% c('CD4', 'CD8')),], aes(x=UMAP1, y=UMAP2), color='lightblue', size=1) +
  geom_point(data=df[which(df$origin == 'donor' & (!df$celltype %in% c('CD4', 'CD8'))),], aes(x=UMAP1, y=UMAP2), color='black', size=1) +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
ggsave('./figures/umaps/20211107_AML1010_donor_derived_cells.png', width = 2.5, height = 2.5)

length(which((!df$celltype %in% c('CD4', 'CD8')) & df$origin != 'none'))
length(which((!df$celltype %in% c('CD4', 'CD8')) & df$origin == 'donor'))

df %>% filter(Sample != 'none') %>% group_by(Sample, celltype) %>% 
  summarize(donor.freq = length(which(origin == 'donor')) / length(origin))


write.csv2(chimerism.df, file = './data/mtDNA/20200302_AML1010_donor_recipient.csv')

