setwd('/Users/shaka87/dfci/asap_seq/')

library(ArchR)
library(ggplot2)
library(Seurat)

CLL.mix.so = readRDS('./data/artifical_mixing/CLL_mix_so.rds')

# read out mixing experiment
df = data.frame()
for (exp in c(1,5,10,50,100,500,1000)) {
  message(exp)
  
  boo = data.table::fread(paste0('./data/artifical_mixing/souporcell/mixing_', exp, '/clusters.tsv'))
  boo$spikein = exp
  df = rbind(df, boo)
}

# doublet rate
ggplot(df %>% group_by(spikein) %>% 
  summarize(doublets = length(which(status == 'doublet'))), aes(x=spikein, y=doublets)) + 
  geom_point(size=0.5) +
  scale_x_log10('CLL1 spiked in') + 
  scale_y_continuous('doublet cells') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/artificial_mixing/20230322_souporcell_doublet_rate.svg', width = 1.5, height = 1.5)

# unassigned rate
ggplot(df %>% group_by(spikein) %>% 
         summarize(unassigned = length(which(status == 'unassigned'))), aes(x=spikein, y=unassigned)) + 
  geom_point(size=0.5) +
  scale_x_log10('CLL1 spiked in') + 
  scale_y_continuous('unassigned cells') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/artificial_mixing/20230322_souporcell_unassignment_rate.svg', width = 1.5, height = 1.5)

# ratio
ggplot(df %>% group_by(spikein) %>% 
         filter(status == 'singlet') %>%
         summarize(ratio = ifelse(length(which(assignment == '0')) < length(which(assignment == '1')),
                                  length(which(assignment == '0')) / length(which(assignment %in% c('0', '1'))),
                                  length(which(assignment == '1')) / length(which(assignment %in% c('0', '1'))))), aes(x=spikein, y=ratio)) + 
  geom_point(size=0.5) +
  scale_x_log10('CLL1 spiked in') + 
  scale_y_continuous('ratio of cell clusters') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/artificial_mixing/20230322_souporcell_cluster_ratio.svg', width = 1.5, height = 1.5)

statistics = data.frame()
for (exp in c(1,5,10,50,100,500,1000)) {
  barcodes1 = data.table::fread(paste0('./data/artifical_mixing/barcodes_scRNAseq/CLL_relapse1_1_barcodes.',exp), header = F)$V1
  barcodes3 = data.table::fread(paste0('./data/artifical_mixing/barcodes_scRNAseq/CLL_relapse3_1_barcodes.',exp), header = F)$V1
  
  boo = df %>% filter(spikein == exp & status == 'singlet')
  statistics = rbind(statistics, 
                     data.frame(spikein = exp,
                                cells.CLL1.0 = length(which(boo$barcode[which(boo$assignment == '0')] %in% barcodes1)),
                                cells.CLL1.1 = length(which(boo$barcode[which(boo$assignment == '1')] %in% barcodes1)),
                                cells.CLL3.0 = length(which(boo$barcode[which(boo$assignment == '0')] %in% barcodes3)),
                                cells.CLL3.1 = length(which(boo$barcode[which(boo$assignment == '1')] %in% barcodes3))))
}
statistics$assignment.1 = apply(statistics, MARGIN = 1, FUN = function(x) {
  if (x['cells.CLL1.0'] > x['cells.CLL1.1']) {
    x['cells.CLL1.0'] / (x['cells.CLL1.0'] + x['cells.CLL1.1'])
  } else {
    x['cells.CLL1.1'] / (x['cells.CLL1.0'] + x['cells.CLL1.1'])
  }
})
statistics$assignment.3 = apply(statistics, MARGIN = 1, FUN = function(x) {
  if (x['cells.CLL3.0'] > x['cells.CLL3.1']) {
    x['cells.CLL3.0'] / (x['cells.CLL3.0'] + x['cells.CLL3.1'])
  } else {
    x['cells.CLL3.1'] / (x['cells.CLL3.0'] + x['cells.CLL3.1'])
  }
})

ggplot() + 
  geom_point(data=statistics, aes(x=spikein, y=100*assignment.1), color='orange', size=0.5) + 
  geom_point(data=statistics, aes(x=spikein, y=100*assignment.3), color='purple', size=0.5) +
  geom_line(data=statistics, aes(x=spikein, y=100*assignment.1), color='orange') + 
  geom_line(data=statistics, aes(x=spikein, y=100*assignment.3), color='purple') +
  scale_x_log10('CLL1 spiked in') + 
  scale_y_continuous('%cluster distribution',limits = c(0,100)) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/artificial_mixing/20230322_souporcell_cluster_distribution.svg', width = 1.5, height = 1.5)