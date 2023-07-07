setwd('/Users/shaka87/dfci/asap_seq/')

library(ArchR)
library(ComplexHeatmap)
library(ggplot2)
library(parallel)
library(Seurat)
library(Signac)

source('./R/variantplot.LP.R')

CLL.mix = loadArchRProject('./data/mixing/singlecell/CLL.mix/')

barcodes.1 = as.data.frame(data.table::fread(file = paste0('./data/mixing/singlecell/barcodes_mtscATACseq/CLL_relapse1_1_barcodes.10'), header = F))
barcodes.3 = as.data.frame(data.table::fread(file = paste0('./data/mixing/singlecell/barcodes_mtscATACseq/CLL_relapse3_1_barcodes.10'), header = F))

informative.variants = as.data.frame(read.csv2(file = './data/mixing/singlecell/informative_variants.csv', sep = '\t'))
rownames(informative.variants) = informative.variants$variant
# load mgatk output
mito.data <- ReadMGATK(dir = paste0('./data/mixing/singlecell/mtscATAC-seq/mixing_10.mgatk/'))

variant.df = IdentifyVariants(mito.data$counts, refallele = mito.data$refallele)

# plot potential maternal variants
high.conf <- variant.df[which(variant.df$n_cells_conf_detected > 0),]
high.conf$pos = 'non.significant'
high.conf$pos[high.conf$strand_correlation > 0.8] = 'significant'

ggplot(data = high.conf, aes(x=strand_correlation, y=log10(vmr), color=pos)) +
  geom_point(size=0.5) +
  labs(x = "Strand concordance", y = "Variance-mean ratio") +
  geom_vline(xintercept = 0.8, color = "black", linetype = 2) +
  scale_color_manual(values = c('non.significant' = "grey", 'significant' = "firebrick", 'germline' = "blue")) +
  scale_y_continuous('log10(VMR)') +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/singlecell/figures/20230518_potential_maternal_variants.svg', width = 1.5, height = 1.5)

variants = variant.df$variant[which(variant.df$strand_correlation > 0.85)]

informative.variants$variant %in% variants

# read out informative mtDNA mutations
VAFs = AlleleFreq(mito.data$counts, variants = variants, assay= NULL)
VAFs.mat = as.data.frame(VAFs)

# perform k means clustering
set.seed(1987)
cluster.variants = kmeans(VAFs.mat, centers = 4)
cluster.cells = kmeans(t(VAFs.mat), centers = 4)

boo = as.data.frame(t(cluster.cells$centers))
#boo = as.data.frame(t(cluster.variants$centers))
colnames(boo) = c('center1', 'center2', 'center3', 'center4')
#boo$individual = ifelse(rownames(boo) %in% barcodes.1$V1, 'CLL1', 'CLL3') 

boo$informative.variant = ifelse(rownames(boo) %in% informative.variants$variant, 'yes', 'no')
boo$individual = informative.variants[rownames(boo), 'individual']
ggplot(boo, aes(x=center3, y=center4, color=individual)) + 
  geom_point(size=0.5) +
  scale_color_manual(values = c('CLL1' = 'orange', 'CLL3' = 'purple')) + 
  theme_classic() + 
  theme(legend.position = 'none') +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/singlecell/figures/20230518_potential_maternal_variants_clustering.svg', width = 1.5, height = 1.5)

# find variant clusters that align best with cell clusters
cluster.df = data.frame()
for (a in seq(1,4)) {
  for (b in seq(1,4)) {
    message(paste(a,b))
    mean = mean(rowMeans(VAFs.mat[names(cluster.variants$cluster)[which(cluster.variants$cluster == a)],
                      names(cluster.cells$cluster)[which(cluster.cells$cluster == b)]]))
    cluster.df = rbind(cluster.df, data.frame(a = a, b = b, mean = mean))
  }
}

# variants 3 and 4 contain maternal variants for clusters 2+4 and 3
informative.variants$variant[which(informative.variants$individual == 'CLL1')] %in% variants.4
informative.variants$variant[which(informative.variants$individual == 'CLL3')] %in% variants.3

variants.1 = names(cluster.variants$cluster)[which(cluster.variants$cluster == 1)]
variants.2 = names(cluster.variants$cluster)[which(cluster.variants$cluster == 2)]
variants.3 = names(cluster.variants$cluster)[which(cluster.variants$cluster == 3)]
variants.4 = names(cluster.variants$cluster)[which(cluster.variants$cluster == 4)]
cells.2 = c(names(cluster.cells$cluster)[which(cluster.cells$cluster == 2)],
            names(cluster.cells$cluster)[which(cluster.cells$cluster == 4)])
cells.3 = names(cluster.cells$cluster)[which(cluster.cells$cluster == 3)]
cells.2 = cells.2[sample(length(cells.2), 1000)]

# clean up variants 2 and 3 to retain only relevant maternal variants
maternal.variants.2 = names(which(rowMeans(VAFs.mat[variants.3, cells.2]) > 0.8 & 
                                    rowMeans(VAFs.mat[variants.3, cells.3]) < 0.2))
maternal.variants.3 = names(which(rowMeans(VAFs.mat[variants.4, cells.3]) > 0.8 & 
                                    rowMeans(VAFs.mat[variants.4, cells.2]) < 0.2))

# sort
maternal.variants.2 = gtools::mixedsort(maternal.variants.2)
maternal.variants.3 = gtools::mixedsort(maternal.variants.3)

# demonstrate result
col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))
svglite::svglite('./figure_mixing/singlecell/figures/20230518_maternal_variants_heatmap.svg', width = 3, height = 3)
Heatmap(VAFs.mat[c(maternal.variants.2, maternal.variants.3),c(cells.2, cells.3)], 
        show_row_dend = F, show_column_dend = F, show_column_names = F, cluster_rows = F,
        col = col_fun, row_names_side = 'left', row_names_gp = gpar(fontsize=6),
        row_split = factor(c(rep('cluster 1', length(maternal.variants.2)),
                             rep('cluster 2', length(maternal.variants.3))),
                           levels = c('cluster 1', 'cluster 2')),
        column_split = factor(c(rep('cluster 1', length(cells.2)),
                                rep('cluster 2', length(cells.3))),
                              levels = c('cluster 1', 'cluster 2')), 
        border = T, use_raster = T, raster_quality = 10)
dev.off()