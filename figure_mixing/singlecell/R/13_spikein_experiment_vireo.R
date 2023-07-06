setwd('/Users/shaka87/dfci/asap_seq/')

library(ArchR)
library(dplyr)
library(ggplot2)
library(Seurat)

CLL.mix.so = readRDS('./data/artifical_mixing/CLL_mix_so.rds')

# read out mixing experiment
df = data.frame()
for (exp in c(1,5,10,50,100,500,1000)) {
  message(exp)
  
  boo = data.table::fread(paste0('./data/artifical_mixing/vireo/mixing_', exp, '/donor_ids.tsv'))
  boo$spikein = exp
  df = rbind(df, boo)
}

# doublet rate
ggplot(df %>% group_by(spikein) %>% 
  summarize(doublets = length(which(donor_id == 'doublet'))), aes(x=spikein, y=doublets)) + 
  geom_point(size=0.5) +
  scale_x_log10('CLL1 spiked in') + 
  scale_y_continuous('doublet cells') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/artificial_mixing/20230327_vireo_doublet_rate.svg', width = 1.5, height = 1.5)

# unassigned rate
ggplot(df %>% group_by(spikein) %>% 
         summarize(unassigned = length(which(donor_id == 'unassigned'))), aes(x=spikein, y=unassigned)) + 
  geom_point(size=0.5) +
  scale_x_log10('CLL1 spiked in') + 
  scale_y_continuous('unassigned cells') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/artificial_mixing/20230327_vireo_unassignment_rate.svg', width = 1.5, height = 1.5)

# identified CLL1 and CLL3 cells
ggplot(data=df %>% group_by(spikein) %>% 
         summarize(CLL1 = length(which(donor_id == 'CLL1')),
                   CLL3 = length(which(donor_id == 'CLL3')))) +
  geom_point(aes(x=spikein, y=CLL1), color='orange', size=0.5) + 
  geom_point(aes(x=spikein, y=CLL3), color='purple', size=0.5) + 
  scale_x_log10('CLL1 spiked in') + 
  scale_y_log10('identified cells') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/artificial_mixing/20230327_vireo_identified_cells.svg', width = 1.5, height = 1.5)

