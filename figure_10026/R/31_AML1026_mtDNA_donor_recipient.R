# perform mtDNA-based donor-recipient deconvolution for AML1026

library(ArchR)
library(gplots)
library(Seurat)
library(ComplexHeatmap)

source('./figure_10026/R/Tcell.colors.R')

AML.1026.mito = loadArchRProject('./data/10026/AML.1026.mito/')
TSB.so = readRDS('./data/10026/objects/20210601_TSB.so.rds')

### 1026
s = 'AML1026'
donor.variants = gtools::mixedsort(c('10463T>C', '16519T>C', '15928G>A', '16126T>C', '16153G>A', '73A>G', '15607A>G', '15452C>A', '7028C>T', '4917A>G',
                                     '13368G>A', '4216T>C', '8697G>A', '11812A>G', '14766C>T', '11719G>A', '709G>A', '8269G>A', '14905G>A', '2706A>G',
                                     '11251A>G', '14233A>G', '1888G>A', '9947G>A', '150C>T', '16294C>T', '16296C>T'))
recipient.variants = gtools::mixedsort(c('16304T>C', '456C>T', '8433T>C', '15833C>T', '4336T>C', '9722T>C', '4011C>T', '93A>G'))


# read mtDNA data
variants = read.table(paste0('./data/10026/mtDNA/20210429_',s,'_high_confidence_variants.csv'), header = T)
combined.frequencies = readRDS(paste0('./data/10026/mtDNA/20210429_',s,'_combined_mutation_frequencies.rds'))

# cells with mtDNA information
cells = AML.1026.mito$cellNames[which(grepl(s,AML.1026.mito$Sample))]
cells = cells[which(cells %in% colnames(combined.frequencies))]
combined.frequencies = combined.frequencies[,cells]

# filter uninformative mtDNA mutations
#combined.frequencies = combined.frequencies[-which(rowSums(combined.frequencies) == 0),]
combined.frequencies = combined.frequencies[-which(rowMeans(combined.frequencies) < 0.01),]
#combined.frequencies = combined.frequencies[variants$variant[which(variants$n_cells_conf_detected > 50)],]
#combined.frequencies = combined.frequencies[order(rowVars(as.matrix(combined.frequencies)), decreasing = T)[1:30],]

# order by phenotype
cells.ordered = c()
for (cluster in names(AML.1026.mito.colors)) {
  boo = AML.1026.mito$cellNames[which(AML.1026.mito$manual.clusters == cluster & AML.1026.mito$cellNames %in% cells)]
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
ggsave(paste0('./figure_10026/figures/deconvolution/20220302_AML1026_mtDNA_variants.svg'), width = 1.5, height = 1.5, plot = p)

df = getEmbedding(AML.1026.mito)
colnames(df) = c('UMAP1', 'UMAP2')
df$origin = chimerism.df[rownames(df),'individual']
df[AML.1026.mito$cellNames,'celltype'] = AML.1026.mito$manual.clusters
df[AML.1026.mito$cellNames,'Sample'] = AML.1026.mito$Sample.hash

ggplot() + 
  geom_point(data=df[which(df$origin == 'recipient'),], aes(x=UMAP1, y=UMAP2), color='lightgrey', size=0.5) +
  geom_point(data=df[which(df$origin == 'donor' & df$celltype %in% c('T cell')),], aes(x=UMAP1, y=UMAP2), color='lightblue', size=1) +
  geom_point(data=df[which(df$origin == 'donor' & (!df$celltype %in% c('T cell'))),], aes(x=UMAP1, y=UMAP2), color='black', size=1) +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
ggsave('./figure_10026/figures/umaps/20211107_AML1026_donor_derived_cells.png', width = 2.5, height = 2.5)

length(which((!df$celltype %in% c('T cell')) & df$origin != 'none'))
length(which((!df$celltype %in% c('T cell')) & df$origin == 'donor'))

df %>% filter(Sample != 'none') %>% group_by(Sample, celltype) %>% 
  summarize(donor.freq = length(which(origin == 'donor')) / length(origin))

write.csv2(chimerism.df, file = './data/10026/mtDNA/20200302_AML1026_donor_recipient.csv')



# score calculated from donor and recipient variants
chimerism.score = as.numeric(colMeans(combined.frequencies[donor.variants,cells.ordered])) - 
  as.numeric(colMeans(combined.frequencies[recipient.variants,cells.ordered]))

chimerism.analysis = ifelse(chimerism.score > 0, 'donor', 'recipient')

plotEmbedding(AML.1026.mito, highlightCells = cells.ordered[which(chimerism.analysis == 'donor')])

df = getEmbedding(AML.1026.mito)
colnames(df) = c('UMAP1', 'UMAP2')
df$origin = 'none'
df[cells.ordered,'origin'] = chimerism.analysis
df[AML.1026.mito$cellNames,'celltype'] = AML.1026.mito$manual.clusters
df[AML.1026.mito$cellNames,'Sample'] = AML.1026.mito$Sample.hash
