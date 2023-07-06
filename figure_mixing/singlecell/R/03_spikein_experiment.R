setwd('/Users/shaka87/dfci/asap_seq/')

library(ArchR)
library(dplyr)
library(ggplot2)
library(parallel)
library(Seurat)
library(Signac)

CLL.mix = loadArchRProject('/Users/shaka87/ArchRProjects/CLL_relapse/CLL.mix/')

informative.variants = as.data.frame(read.csv2(file = './data/artifical_mixing/informative_variants.csv', sep = '\t'))

# read out mixing experiment
df = data.frame()
for (exp in c(1,5,10,50,100,500,1000)) {
  message(exp)
  # load mgatk output
  mito.data <- ReadMGATK(dir = paste0('./data/artifical_mixing/mtscATAC-seq/mixing_',exp,'.mgatk/'))
  
  # read out informative mtDNA mutations
  VAFs = AlleleFreq(mito.data$counts, variants = informative.variants$variant, assay= NULL)
  
  df = rbind(df, data.frame(bc = colnames(VAFs),
                            VAFs.1 = colMeans(VAFs[informative.variants$variant[which(informative.variants$individual == 'CLL1')],]),
                            VAFs.3 = colMeans(VAFs[informative.variants$variant[which(informative.variants$individual == 'CLL3')],]),
                            spike.in = exp))
}
write.table(df, quote = F, sep = '\t', row.names = F, file = './data/artifical_mixing/mtscATAC-seq/results.csv')
df = read.table(file = './data/artifical_mixing/mtscATAC-seq/results.csv', header = T)

# read ground-truth data
df$sample.orig = 'none'
for (exp in c(1,5,10,50,100,500,1000)) {
  barcodes.1 = as.data.frame(data.table::fread(file = paste0('./data/artifical_mixing/barcodes_mtscATACseq/CLL_relapse1_1_barcodes.', exp), header = F))
  barcodes.3 = as.data.frame(data.table::fread(file = paste0('./data/artifical_mixing/barcodes_mtscATACseq/CLL_relapse3_1_barcodes.', exp), header = F))
  df[which(df$spike.in == exp & df$bc %in% barcodes.1$V1), 'sample.orig'] = 'CLL1'
  df[which(df$spike.in == exp & df$bc %in% barcodes.3$V1), 'sample.orig'] = 'CLL3'
}

# annotate based on mtDNA
df$sample.annotation = 'none'
df[which(df$VAFs.1 > 0.8 & df$VAFs.3 < 0.2), 'sample.annotation'] = 'CLL1'
df[which(df$VAFs.1 < 0.8 & df$VAFs.3 > 0.8), 'sample.annotation'] = 'CLL3'

# assess annotation
df$result = apply(df, MARGIN = 1, FUN = function(x) {
  if (x['sample.annotation'] == 'none') {
    'none'
  } else if (x['sample.orig'] == x['sample.annotation']) {
    'correct'
  } else if (x['sample.orig'] != x['sample.annotation']) {
    'mismatch'
}})

results = df %>% group_by(spike.in) %>% 
  summarize(none.1 = length(which(sample.orig == 'CLL1' & result == 'none')),
            none.3 = length(which(sample.orig == 'CLL3' & result == 'none')),
            correct.1 = length(which(sample.orig == 'CLL1' & result == 'correct')),
            correct.3 = length(which(sample.orig == 'CLL3' & result == 'correct'))) %>% 
  as.data.frame()

ggplot() + 
  geom_abline(slope = 1) + 
  geom_point(data=results, aes(x=spike.in, y=correct.1), size=0.5, color='orange') + 
  geom_point(data=results, aes(x=spike.in, y=correct.3), size=0.5, color='purple') + 
  scale_x_log10('cells spiked-in') + 
  scale_y_log10('correct annotation') +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/artificial_mixing/20230322_spikein_mtscATAC.svg', width = 1.5, height = 1.5)

ggplot(results, aes(x=spike.in, y=none.1)) + 
  geom_abline(slope = 1) + 
  geom_point(size=0.5) + 
  scale_x_log10('cells spiked-in', limits = c(1,1000)) + 
  scale_y_log10('missing annotation', limits = c(1,1000)) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/artificial_mixing/20230322_spikein_mtscATAC_unannotated.svg', width = 1.5, height = 1.5)

