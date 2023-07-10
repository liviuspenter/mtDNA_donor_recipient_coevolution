# create ArchR object just for patient IST1 (IST1_1, IST1_2, IST2_1, IST2_2) and explore mtDNA mutations

library(ArchR)
library(ComplexHeatmap)
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
combined.mutation.frequencies = readRDS('./data/IST/mtDNA/20220110_IST1_2_combined_mutation_frequencies.rds')
combined.mutation.frequencies.12 = readRDS(file = './data/IST/mtDNA/20220110_IST1_2_combined_mutation_frequencies.rds')

germline.variants = read.csv2(file='./data/IST/objects/20220117_IST_germline_variants.csv')
# germline variants with incomplete coverage
exclude.variants = c('310T>C', '5458T>C', '7457G>A')
non.germline.variants = setdiff(rownames(combined.mutation.frequencies), c(exclude.variants, germline.variants$variant[which(germline.variants$sample == 'IST1')]))

IST.asap.mito.1 = subsetArchRProject(ArchRProj = IST.asap.mito, 
                                     cells = IST.asap.mito$cellNames[which(IST.asap.mito$patient == 'IST1')],
                                     outputDirectory = './data/IST/IST.asap.mito.1/', threads = 12, force = T)
IST.asap.mito.1 = addClusters(input = IST.asap.mito.1, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 0.3, force = T)
IST.asap.mito.1 = addUMAP(ArchRProj = IST.asap.mito.1, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
IST.asap.mito.1 = addImputeWeights(IST.asap.mito.1)

IST.asap.mito.1 = loadArchRProject('./data/IST/IST.asap.mito.1/')

# find recipient-specific mtDNA mutations
# 3919T>C and 10776T>C
mtDNA.so.recipient = CreateSeuratObject(combined.mutation.frequencies[non.germline.variants,IST.asap.mito.1$cellNames[which(IST.asap.mito.1$individual == 'recipient')]])
mtDNA.so.recipient = FindVariableFeatures(mtDNA.so.recipient)
mtDNA.so.recipient = ScaleData(mtDNA.so.recipient, features = rownames(mtDNA.so.recipient))
mtDNA.so.recipient = RunPCA(mtDNA.so.recipient)
mtDNA.so.recipient = RunUMAP(mtDNA.so.recipient, dims = 1:30)
mtDNA.so.recipient = FindNeighbors(mtDNA.so.recipient)
mtDNA.so.recipient = FindClusters(mtDNA.so.recipient, resolution = 0.2)
mtDNA.so.recipient$Sample = stringr::str_split_fixed(colnames(mtDNA.so.recipient), pattern = '#', n=2)[,1]

mtDNA.markers <- FindAllMarkers(mtDNA.so.recipient, min.pct = 0.05, logfc.threshold = 0.01)
mtDNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(mtDNA.so.recipient), group.by = c('seurat_clusters', 'Sample'), features = top10$gene, slot = 'counts', disp.max = 0.1, 
          group.colors =  BuenColors::jdb_palette(name = 'corona'), label = F, raster = T) + #NoLegend() + 
  scale_fill_gradientn(colours = c('white', BuenColors::jdb_palette(name = 'solar_rojos')), na.value = 'black') +
  theme(axis.text = element_text('Arial', size=8, color='black'))

# find donor-specific mtDNA mutations
# 13326T>C, 8995G>A, 1793G>A, 2623A>G, 13785C>T
mtDNA.so.donor = CreateSeuratObject(combined.mutation.frequencies[non.germline.variants,IST.asap.mito.1$cellNames[which(IST.asap.mito.1$individual == 'donor')]])
mtDNA.so.donor = FindVariableFeatures(mtDNA.so.donor)
mtDNA.so.donor = ScaleData(mtDNA.so.donor, features = rownames(mtDNA.so.donor))
mtDNA.so.donor = RunPCA(mtDNA.so.donor)
mtDNA.so.donor = RunUMAP(mtDNA.so.donor, dims = 1:30)
mtDNA.so.donor = FindNeighbors(mtDNA.so.donor)
mtDNA.so.donor = FindClusters(mtDNA.so.donor, resolution = 0.2)
mtDNA.so.donor$Sample = stringr::str_split_fixed(colnames(mtDNA.so.donor), pattern = '#', n=2)[,1]

mtDNA.markers <- FindAllMarkers(mtDNA.so.donor, min.pct = 0.05, logfc.threshold = 0.01)
mtDNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(mtDNA.so.donor), group.by = c('seurat_clusters', 'Sample'), features = top10$gene, slot = 'counts', disp.max = 0.1, 
          group.colors =  BuenColors::jdb_palette(name = 'corona'), label = F, raster = T) + #NoLegend() + 
  scale_fill_gradientn(colours = c('white', BuenColors::jdb_palette(name = 'solar_rojos')), na.value = 'black') +
  theme(axis.text = element_text('Arial', size=8, color='black'))

# plot individual mtDNA mutations
df = data.frame(bc = IST.asap.mito.1$cellNames,
                sample = IST.asap.mito.1$Sample,
                individual = IST.asap.mito.1$individual,
                manual.cluster = IST.asap.mito.1$manual.cluster,
                UMAP1 = getEmbedding(IST.asap.mito.1)[,1],
                UMAP2 = getEmbedding(IST.asap.mito.1)[,2])

df = cbind(df, t(combined.mutation.frequencies.12[c('3919T>C', '5458T>C', '7457G>A', '10776T>C', '1793G>A', '2623A>G', '8995G>A','13326T>C','13785C>T'), df$bc]))

df.long = df %>%
  tidyr::pivot_longer(cols = c('3919T>C', '10776T>C', '13326T>C', '8995G>A', '1793G>A', '2623A>G', '13785C>T'),
                      names_to = 'variant', values_to = 'heteroplasmy') 
  
boo = df.long %>% group_by(sample, individual, manual.cluster, variant) %>% summarize(detected = length(which(heteroplasmy != 0)),
                                                                                      mean = mean(heteroplasmy),
                                                                                      n = length(manual.cluster),
                                                                                      freq = length(which(heteroplasmy != 0)) / 
                                                                                        length(manual.cluster))

col_fun = circlize::colorRamp2(breaks = seq(0,0.1,0.1/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))

cells = df$bc[which(df$individual == 'donor' & df$sample == 'IST1_1')]
ha = columnAnnotation(celltype = factor(df[cells, 'manual.cluster'], levels = c('HSC', 'Myeloid', 'Erythroid','TNK', 'B cell')))
Heatmap(combined.mutation.frequencies.12[c('3919T>C', '10776T>C', '13326T>C', '8995G>A', '1793G>A', '2623A>G', '13785C>T'),cells],
        show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F,
        col = col_fun, column_split = df[cells, 'manual.cluster'], border=T, top_annotation = ha)

cells = df$bc[which(df$individual == 'donor' & df$sample == 'IST1_2')]
ha = columnAnnotation(celltype = factor(df[cells, 'manual.cluster'], levels = c('HSC', 'Myeloid', 'Erythroid','TNK', 'B cell')))
Heatmap(combined.mutation.frequencies.12[c('3919T>C', '10776T>C', '13326T>C', '8995G>A', '1793G>A', '2623A>G', '13785C>T'),cells],
        show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F,
        col = col_fun, column_split = df[cells, 'manual.cluster'], border=T, top_annotation = ha)

cells = df$bc[which(df$individual == 'donor' & df$sample == 'IST2_1')]
ha = columnAnnotation(celltype = factor(df[cells, 'manual.cluster'], levels = c('HSC', 'Myeloid', 'Erythroid','TNK', 'B cell')))
Heatmap(combined.mutation.frequencies.12[c('3919T>C', '10776T>C', '13326T>C', '8995G>A', '1793G>A', '2623A>G', '13785C>T'),cells],
        show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F,
        col = col_fun, column_split = df[cells, 'manual.cluster'], border=T, top_annotation = ha)

cells = df$bc[which(df$individual == 'donor' & df$sample == 'IST2_2')]
ha = columnAnnotation(celltype = factor(df[cells, 'manual.cluster'], levels = c('HSC', 'Myeloid', 'Erythroid','TNK', 'B cell')))
Heatmap(combined.mutation.frequencies.12[c('3919T>C', '10776T>C', '13326T>C', '8995G>A', '1793G>A', '2623A>G', '13785C>T'),cells],
        show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F,
        col = col_fun, column_split = df[cells, 'manual.cluster'], border=T, top_annotation = ha)

for (variant in c('3919T>C', '5458T>C', '7457G>A', '10776T>C', '1793G>A', '2623A>G', '8995G>A','13326T>C','13785C>T')) {
  p=ggplot(df[order(df[,variant]),], aes(x=UMAP1, y=UMAP2)) + geom_point(size=1, aes(color=df[order(df[,variant]),variant])) + 
    scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'brewer_orange'))) + 
    theme_classic() +
    theme(legend.position = 'none') +
    NoAxes()
  ggsave(paste0('./figure_IST/figures/IST1/UMAPs/20230404_IST1_', variant, '.png'), width = 3, height = 3, dpi = 600, plot = p)
}

# plot mitochondrial barcodes
for (mutation in rownames(combined.mutation.frequencies)) {
  IST.asap.mito.1$vaf = unlist(combined.mutation.frequencies[mutation, IST.asap.mito.1$cellNames])
  IST.asap.mito.1$vaf[which(IST.asap.mito.1$vaf > 0.1)] = 0.1
  p= plotEmbedding(IST.asap.mito.1, name = 'vaf', pal = c('grey','red'), plotAs = 'points', na.rm=T) + 
    ggtitle(mutation) + theme(plot.title = element_text(hjust = 0.5))
  plotPDF(p, name = paste0('IST1.mito.',mutation), ArchRProj = IST.asap.mito.1, addDOC = F, width = 4, height = 4)
}

# % cells marked by mtDNA mutations in donor and recipient, pre and post per celltype
donor.cells = IST.asap.mito.1$cellNames[which(IST.asap.mito.1$individual == 'donor')]
recipient.cells = IST.asap.mito.1$cellNames[which(IST.asap.mito.1$individual == 'recipient')]
mtDNA.statistics.all = data.frame()
mtDNA.statistics = list()
for (celltype in unique(IST.asap.mito.1$manual.cluster)) {
  mtDNA.statistics = append(mtDNA.statistics, list(data.frame(mtDNA.mutation = as.character(), sample = as.character(), 
                                                              cells = as.numeric(), frequency = as.numeric())))
}
names(mtDNA.statistics) = unique(IST.asap.mito.1$manual.cluster)

# number of mtDNA mutations per celltype
for (sample in unique(IST.asap.mito.1$Sample)) {
  for (celltype in unique(IST.asap.mito.1$manual.cluster)) {
    cells = intersect(IST.asap.mito.1$cellNames[which(IST.asap.mito.1$Sample == sample & IST.asap.mito.1$manual.cluster == celltype)], donor.cells)
    cells.2 = intersect(IST.asap.mito.1$cellNames[which(IST.asap.mito.1$Sample == sample & IST.asap.mito.1$manual.cluster == celltype)], recipient.cells)
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
write.csv2(mtDNA.statistics.all, file = './data/IST/mtDNA/20220119_IST1_mtDNA_statistics.csv', quote = F, row.names = F)
