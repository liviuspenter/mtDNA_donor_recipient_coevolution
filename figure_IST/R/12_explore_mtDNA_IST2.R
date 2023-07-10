# create ArchR object just for patient IST2 (samples IST3_1, IST3_2) and explore mtDNA mutations

library(ArchR)
library(ggplot2)
library(ggrepel)
library(gplots)
library(grid)
library(parallel)
library(Signac)
library(Seurat)
library(BuenColors)
library(dplyr)

IST.asap.mito = loadArchRProject('./data/IST/IST.asap.mito/')
combined.mutation.frequencies = readRDS('./data/IST/mtDNA/20220110_IST3_combined_mutation_frequencies.rds')

germline.variants = read.csv2(file='./data/IST/objects/20220117_IST_germline_variants.csv')
# germline variants with incomplete coverage
exclude.variants = c('310T>C')
non.germline.variants = setdiff(rownames(combined.mutation.frequencies), c(exclude.variants, germline.variants$variant[which(germline.variants$sample == 'IST3')]))

IST.asap.mito.3 = subsetArchRProject(ArchRProj = IST.asap.mito, 
                                     cells = IST.asap.mito$cellNames[which(IST.asap.mito$patient == 'IST3')],
                                     outputDirectory = './data/IST/IST.asap.mito.3/', threads = 12, force = T)
IST.asap.mito.3 = addClusters(input = IST.asap.mito.3, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 0.3, force = T)
IST.asap.mito.3 = addUMAP(ArchRProj = IST.asap.mito.3, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
IST.asap.mito.3 = addHarmony(IST.asap.mito.3)
IST.asap.mito.3$manual.cluster[which(IST.asap.mito.3$Clusters %in% c('C1'))] = 'Erythroid'
IST.asap.mito.3$manual.cluster[which(IST.asap.mito.3$Clusters %in% c('C11', 'C12', 'C13'))] = 'TNK'
IST.asap.mito.3$manual.cluster[which(IST.asap.mito.3$Clusters %in% c('C10'))] = 'HSC'
IST.asap.mito.3$manual.cluster[which(IST.asap.mito.3$Clusters %in% c('C8', 'C9'))] = 'B cell'
IST.asap.mito.3$manual.cluster[which(IST.asap.mito.3$Clusters %in% c('C2','C3', 'C4', 'C5', 'C6', 'C7'))] = 'Myeloid'
IST.asap.mito.3 = addImputeWeights(IST.asap.mito.3)
saveArchRProject(IST.asap.mito.3)

IST.asap.mito.3 = loadArchRProject('./data/IST/IST.asap.mito.3/')


# find donor-specific mtDNA mutations
# 
mtDNA.so.donor = CreateSeuratObject(combined.mutation.frequencies[non.germline.variants,IST.asap.mito.3$cellNames[which(IST.asap.mito.3$individual == 'donor')]])
mtDNA.so.donor = FindVariableFeatures(mtDNA.so.donor)
mtDNA.so.donor = ScaleData(mtDNA.so.donor, features = rownames(mtDNA.so.donor))
mtDNA.so.donor = RunPCA(mtDNA.so.donor)
mtDNA.so.donor = RunUMAP(mtDNA.so.donor, dims = 1:30)
mtDNA.so.donor = FindNeighbors(mtDNA.so.donor)
mtDNA.so.donor = FindClusters(mtDNA.so.donor, resolution = 0.3)
mtDNA.so.donor$Sample = stringr::str_split_fixed(colnames(mtDNA.so.donor), pattern = '#', n=2)[,1]

mtDNA.markers <- FindAllMarkers(mtDNA.so.donor, min.pct = 0.05, logfc.threshold = 0.01)
mtDNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(mtDNA.so.donor), group.by = c('seurat_clusters', 'Sample'), features = top10$gene, slot = 'counts', disp.max = 0.1, 
          group.colors =  BuenColors::jdb_palette(name = 'corona'), label = F, raster = T) + #NoLegend() + 
  scale_fill_gradientn(colours = c('white', BuenColors::jdb_palette(name = 'solar_rojos')), na.value = 'black') +
  theme(axis.text = element_text('Arial', size=8, color='black'))


# plot mitochondrial barcodes
for (mutation in rownames(combined.mutation.frequencies)) {
  IST.asap.mito.3$vaf = unlist(combined.mutation.frequencies[mutation, IST.asap.mito.3$cellNames])
  IST.asap.mito.3$vaf[which(IST.asap.mito.3$vaf > 0.1)] = 0.1
  p= plotEmbedding(IST.asap.mito.3, name = 'vaf', pal = c('grey','red'), plotAs = 'points', na.rm=T) + 
    ggtitle(mutation) + theme(plot.title = element_text(hjust = 0.5))
  plotPDF(p, name = paste0('IST3.mito.',mutation), ArchRProj = IST.asap.mito.3, addDOC = F, width = 4, height = 4)
}

# % cells marked by mtDNA mutations in donor and recipient, pre and post per celltype
donor.cells = IST.asap.mito.3$cellNames[which(IST.asap.mito.3$individual == 'donor')]
recipient.cells = IST.asap.mito.3$cellNames[which(IST.asap.mito.3$individual == 'recipient')]
mtDNA.statistics.all = data.frame()
mtDNA.statistics = list()
for (celltype in unique(IST.asap.mito.3$manual.cluster)) {
  mtDNA.statistics = append(mtDNA.statistics, list(data.frame(mtDNA.mutation = as.character(), sample = as.character(), 
                                                              cells = as.numeric(), frequency = as.numeric())))
}
names(mtDNA.statistics) = unique(IST.asap.mito.3$manual.cluster)

# number of mtDNA mutations per celltype
for (sample in unique(IST.asap.mito.3$Sample)) {
  for (celltype in unique(IST.asap.mito.3$manual.cluster)) {
    cells = intersect(IST.asap.mito.3$cellNames[which(IST.asap.mito.3$Sample == sample & IST.asap.mito.3$manual.cluster == celltype)], donor.cells)
    cells.2 = intersect(IST.asap.mito.3$cellNames[which(IST.asap.mito.3$Sample == sample & IST.asap.mito.3$manual.cluster == celltype)], recipient.cells)
    for (mtDNA.mutation in non.germline.variants) {
      mtDNA.statistics[celltype][[1]] = rbind(mtDNA.statistics[celltype][[1]], 
                                         data.frame(mtDNA.mutation = mtDNA.mutation, sample = sample, 
                                                    cells = length(which(combined.mutation.frequencies[mtDNA.mutation,cells] > 0.05)),
                                                    frequency = length(which(combined.mutation.frequencies[mtDNA.mutation,cells] > 0.05)) / length(cells)))
      mtDNA.statistics.all = rbind(mtDNA.statistics.all, data.frame(mtDNA.mutation = mtDNA.mutation, 
                                                                    sample = sample, 
                                                                    cells = length(which(combined.mutation.frequencies[mtDNA.mutation,cells] > 0.05)),
                                                                    frequency = length(which(combined.mutation.frequencies[mtDNA.mutation,cells] > 0.05)) / 
                                                                      length(cells),
                                                                    celltype = celltype, 
                                                                    individual = 'donor'))
      mtDNA.statistics.all = rbind(mtDNA.statistics.all, data.frame(mtDNA.mutation = mtDNA.mutation, 
                                                                    sample = sample, 
                                                                    cells = length(which(combined.mutation.frequencies[mtDNA.mutation,cells.2] > 0.05)),
                                                                    frequency = length(which(combined.mutation.frequencies[mtDNA.mutation,cells.2] > 0.05)) / 
                                                                      length(cells.2),
                                                                    celltype = celltype, 
                                                                    individual = 'recipient'))
    }
  }
}
write.csv2(mtDNA.statistics.all, file = './data/IST/mtDNA/20220119_IST3_mtDNA_statistics.csv', quote = F, row.names = F)

# plot individual mtDNA mutations
df = data.frame(bc = IST.asap.mito.3$cellNames,
                sample = IST.asap.mito.3$Sample,
                individual = IST.asap.mito.3$individual,
                manual.cluster = IST.asap.mito.3$manual.cluster,
                UMAP1 = getEmbedding(IST.asap.mito.3)[,1],
                UMAP2 = getEmbedding(IST.asap.mito.3)[,2])

df$UMAP1 = getEmbedding(IST.asap.mito)[df$bc,1]
df$UMAP2 = getEmbedding(IST.asap.mito)[df$bc,2]

df = cbind(df, t(combined.mutation.frequencies[c('6876G>A', '6701A>G','10290G>C'), df$bc]))

for (variant in c('6876G>A', '6701A>G', '10290G>C')) {
  p=ggplot(df[order(df[,variant]),], aes(x=UMAP1, y=UMAP2)) + 
    geom_point(size=1, aes(color=df[order(df[,variant]),variant])) + 
    scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'brewer_orange'))) + 
    theme_classic() +
    theme(legend.position = 'none') +
    NoAxes()
  ggsave(paste0('./figure_IST/figures/IST2/UMAPs/20230626_IST2_', variant, '.png'), width = 3, height = 3, dpi = 600, plot = p)
}