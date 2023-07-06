setwd('/Users/shaka87/dfci/asap_seq/')

library(ArchR)
library(ggplot2)
library(parallel)
library(Seurat)
library(Signac)

source('./analysis/variantplot.LP.R')

CLL.mix = loadArchRProject('/Users/shaka87/ArchRProjects/CLL_relapse/CLL.mix/')

# extract mtDNA variants from "germline" controls
mito.data.1 = ReadMGATK('./data/artifical_mixing/CLL_relapse1_1.mgatk/final/')
mito.data.3 = ReadMGATK('./data/artifical_mixing/CLL_relapse3_1.mgatk/final/')
variant.df.1 = IdentifyVariants(mito.data.1$counts, refallele = mito.data.1$refallele)
write.table(variant.df.1, quote = F, sep = '\t', row.names = F, file = './data/artifical_mixing/CLL_relapse1_1.variants.csv')
variant.df.3 = IdentifyVariants(mito.data.3$counts, refallele = mito.data.3$refallele)
write.table(variant.df.3, quote = F, sep = '\t', row.names = F, file = './data/artifical_mixing/CLL_relapse3_1.variants.csv')

VariantPlot.LP(variant.df.1)
ggsave('./figures/artificial_mixing/20230321_germline1.svg', width = 1.5, height = 1.5)
VariantPlot.LP(variant.df.3)
ggsave('./figures/artificial_mixing/20230321_germline3.svg', width = 1.5, height = 1.5)

variants.1 = variant.df.1$variant[which(variant.df.1$vmr < 0.01 & variant.df.1$strand_correlation > 0.65)]
variants.3 = variant.df.3$variant[which(variant.df.3$vmr < 0.01 & variant.df.3$strand_correlation > 0.65)]
informative.variants.1 = setdiff(variants.1, variants.3)
informative.variants.3 = setdiff(variants.3, variants.1)

informative.variants = data.frame(variant = c(informative.variants.1, informative.variants.3))
informative.variants$individual = ifelse(informative.variants$variant %in% variants.1, 'CLL1', 'CLL3')

write.table(informative.variants, quote = F, sep = '\t', row.names = F, file = './data/artifical_mixing/informative_variants.csv')

# test for false-positive calls in opposite germline and unannotated cells
VAFs.1 = AlleleFreq(mito.data.1$counts, variants = informative.variants$variant, assay= NULL)
colnames(VAFs.1) = paste0('CLL_relapse1_1#', colnames(VAFs.1))
VAFs.3 = AlleleFreq(mito.data.3$counts, variants = informative.variants$variant, assay= NULL)
colnames(VAFs.3) = paste0('CLL_relapse3_1#', colnames(VAFs.3))

df = data.frame(bc = c(colnames(VAFs.1), colnames(VAFs.3)),
                variants.1 = c(colMeans(VAFs.1[informative.variants.1, ]), 
                               colMeans(VAFs.3[informative.variants.1, ])),
                variants.3 = c(colMeans(VAFs.1[informative.variants.3, ]), 
                               colMeans(VAFs.3[informative.variants.3, ])),
                sample = c(rep('CLL1', ncol(VAFs.1)), rep('CLL3', ncol(VAFs.3))))

df = df[intersect(df$bc, rownames(CLL.mix)),]

df$coverage = c(mito.data.1$depth[stringr::str_split_fixed(df$bc[which(df$sample == 'CLL1')], pattern = '#', n=2)[,2],],
                mito.data.3$depth[stringr::str_split_fixed(df$bc[which(df$sample == 'CLL3')], pattern = '#', n=2)[,2],])
df$individual = 'none'
df[which(df$variants.1 > 0.8 & df$variants.3 < 0.2), 'individual'] = 'CLL1'
df[which(df$variants.3 > 0.8 & df$variants.1 < 0.2), 'individual'] = 'CLL3'

# get information for UMAPs from ArchR object
df = cbind(df, getEmbedding(CLL.mix)[df$bc,])
df$manual.cluster = CLL.mix[df$bc]$manual.cluster

# plot stuff

ggplot(df, aes(x=100*variants.1, y=100*variants.3)) + 
  geom_hline(yintercept = c(20, 80)) + 
  geom_vline(xintercept = c(20, 80)) + 
  ggrastr::rasterize(geom_point(aes(color=individual), size=0.5), dpi=600) +
  scale_x_continuous('% informative variants 1') +
  scale_y_continuous('% informative variants 3') +
  scale_color_manual(values = c('none' = 'black', 'CLL1' = 'orange', 'CLL3' = 'purple')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/artificial_mixing/20230321_VAF_plot.svg', width = 2, height = 2)

p=ggplot(df, aes(x=individual, y=coverage)) + 
  ggrastr::rasterize(geom_jitter(size=0.5, color='grey'), dpi=600) + 
  geom_boxplot(outlier.color = NA, fill=NA, aes(color=individual)) +
  scale_x_discrete(labels = c('CLL1', 'CLL3', 'unannotated')) + 
  scale_y_log10('mtDNA coverage') +
  scale_color_manual(values = c('none' = 'black', 'CLL1' = 'orange', 'CLL3' = 'purple')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
ggsave('./figures/artificial_mixing/20230321_mtDNA_coverage.svg', width = 1.5, height = 2, plot = p)

ggplot() +
  geom_point(data=df, size=0.5,
             aes(x=`IterativeLSI#UMAP_Dimension_1`, y=`IterativeLSI#UMAP_Dimension_2`, color=manual.cluster)) + 
  scale_color_manual(values = c('Mono' = 'darkgreen', 'T cell' = 'blue', 'CLL1' = 'orange', 'CLL3' = 'purple')) +
  theme_classic() + 
  theme(legend.position = 'none') + NoAxes()
ggsave('./figures/artificial_mixing/20230321_UMAP.png', width = 4, height = 4, dpi = 600)

ggplot() +
  geom_point(data=df, size=0.5,
             aes(x=`IterativeLSI#UMAP_Dimension_1`, y=`IterativeLSI#UMAP_Dimension_2`, color=manual.cluster), alpha=0.05) + 
  geom_point(data=df[which(df$individual == 'none'),], size=1,
             aes(x=`IterativeLSI#UMAP_Dimension_1`, y=`IterativeLSI#UMAP_Dimension_2`), color='black') +
  scale_color_manual(values = c('Mono' = 'darkgreen', 'T cell' = 'blue', 'CLL1' = 'orange', 'CLL3' = 'purple')) +
  theme_classic() + 
  theme(legend.position = 'none') + NoAxes()
ggsave('./figures/artificial_mixing/20230321_UMAP_not_annotated.png', width = 4, height = 4, dpi = 600)
