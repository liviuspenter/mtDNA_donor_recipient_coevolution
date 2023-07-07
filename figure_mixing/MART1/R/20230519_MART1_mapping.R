setwd('/Users/liviuspenter/dfci/asap_seq/')

library(ArchR)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(parallel)
library(Seurat)

TSAB.so = readRDS(file='./data/objects/20220225_MART1_TSAB.so')

MART1.mito = loadArchRProject('/Users/shaka87/ArchRProjects/MART1/MART1.mito/')

Pool108.1.bc = data.table::fread('./data/MART1_1/filtered_peak_bc_matrix/barcodes.tsv', header = F)
Pool108.2.bc = data.table::fread('./data/MART1_2/filtered_peak_bc_matrix/barcodes.tsv', header = F)
Pool108.3.bc = data.table::fread('./data/MART1_3/filtered_peak_bc_matrix/barcodes.tsv', header = F)
Pool108.1.bc$V1 = paste0('MART1_1#', Pool108.1.bc$V1)
Pool108.2.bc$V1 = paste0('MART1_1#', Pool108.2.bc$V1)
Pool108.3.bc$V1 = paste0('MART1_1#', Pool108.3.bc$V1)

Pool108.1 = rbind(data.table::fread('./data/MART1_1/Pool108_1_S1_L001.bed'),
                  data.table::fread('./data/MART1_1/Pool108_1_S1_L002.bed'))
colnames(Pool108.1) = c('gene', 'start', 'end', 'bc', 'n')
Pool108.1$bc = as.character(reverseComplement(DNAStringSet(Pool108.1$bc)))
Pool108.1$sample = 'Pool108_1'
Pool108.1$bc = paste0('MART1_1#', Pool108.1$bc, '-1')

Pool108.2 = rbind(data.table::fread('./data/MART1_2/Pool108_2_S2_L001.bed'),
                  data.table::fread('./data/MART1_2/Pool108_2_S2_L002.bed'))
colnames(Pool108.2) = c('gene', 'start', 'end', 'bc', 'n')
Pool108.2$bc = as.character(reverseComplement(DNAStringSet(Pool108.2$bc)))
Pool108.2$sample = 'Pool108_2'
Pool108.2$bc = paste0('MART1_2#', Pool108.2$bc, '-1')

Pool108.3 = rbind(data.table::fread('./data/MART1_3/Pool108_3_S3_L001.bed'),
                  data.table::fread('./data/MART1_3/Pool108_3_S3_L002.bed'))
colnames(Pool108.3) = c('gene', 'start', 'end', 'bc', 'n')
Pool108.3$bc = as.character(reverseComplement(DNAStringSet(Pool108.3$bc)))
Pool108.3$sample = 'Pool108_3'
Pool108.3$bc = paste0('MART1_3#', Pool108.3$bc, '-1')

TCR.df = bind_rows(Pool108.1, Pool108.2, Pool108.3) %>% 
  filter(gene != 'EF1a_promoter') %>%
  filter(bc %in% MART1.mito$cellNames)

df = data.frame(bc = MART1.mito$cellNames,
                individual = MART1.mito$donor)
rownames(df) = df$bc

TCR.df$individual = df[TCR.df$bc, 'individual']

df = cbind(df, getEmbedding(MART1.mito))
colnames(df) = c('bc', 'individual', 'UMAP1', 'UMAP2')

ggplot() + 
  geom_point(data=df, aes(x=UMAP1, y=UMAP2), size=0.5, color='grey90') + 
  geom_point(data=df[which(df$bc %in% TCR.df$bc),], aes(x=UMAP1, y=UMAP2), size=1, color='firebrick') + 
  theme_classic() + 
  NoAxes()
ggsave('./figures/MART1/20230519_MART1_TCR.png', width = 4, height = 4, dpi = 600)