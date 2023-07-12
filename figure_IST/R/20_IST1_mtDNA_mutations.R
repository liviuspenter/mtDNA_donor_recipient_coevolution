library(dplyr)
library(ggplot2)

### 
mtDNA.statistics.all.1 = read.csv2(file = './data/IST/mtDNA/20220119_IST1_mtDNA_statistics.csv')

mtDNA.statistics.all.1 = reshape2::dcast(data = mtDNA.statistics.all.1, 
                                         formula = mtDNA.mutation ~ celltype + sample + individual, 
                                         value.var = 'frequency')

boo = mtDNA.statistics.all.1[, c('mtDNA.mutation', 'TNK_IST2_1_donor', 'TNK_IST2_2_donor')] %>%
  tidyr::pivot_longer(cols = c('TNK_IST2_1_donor', 'TNK_IST2_2_donor'), names_to = 'timepoint')
boo$value[which(boo$value == 0)] = 0.001
ggplot(boo, aes(x=timepoint, y=100*value, group=mtDNA.mutation)) + 
  geom_line(color='grey') +
  geom_point(size=0.5) + 
  scale_x_discrete(labels = c('pre-IST', 'post-IST')) + 
  scale_y_log10('% T/NK cells', limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100), labels = c('ND', '1', '10','100')) +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
ggsave('./figure_IST/figures/IST1/plots/20230711_mtDNA_dynamics_BM_TNK.svg', width = 1.1, height = 2)

boo = mtDNA.statistics.all.1[, c('mtDNA.mutation', 'B cell_IST2_1_donor', 'B cell_IST2_2_donor')] %>%
  tidyr::pivot_longer(cols = c('B cell_IST2_1_donor', 'B cell_IST2_2_donor'), names_to = 'timepoint')
boo$value[which(boo$value == 0)] = 0.001
ggplot(boo, aes(x=timepoint, y=100*value, group=mtDNA.mutation)) + 
  geom_line(color='grey') +
  geom_point(size=0.5) + 
  scale_x_discrete(labels = c('pre-IST', 'post-IST')) + 
  scale_y_log10('% T/NK cells', limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100), labels = c('ND', '1', '10','100')) +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
ggsave('./figure_IST/figures/IST1/plots/20230711_mtDNA_dynamics_BM_Bcell.svg', width = 1.1, height = 2)

boo = mtDNA.statistics.all.1[, c('mtDNA.mutation', 'Myeloid_IST2_1_donor', 'Myeloid_IST2_2_donor')] %>%
  tidyr::pivot_longer(cols = c('Myeloid_IST2_1_donor', 'Myeloid_IST2_2_donor'), names_to = 'timepoint')
boo$value[which(boo$value == 0)] = 0.001
ggplot(boo, aes(x=timepoint, y=100*value, group=mtDNA.mutation)) + 
  geom_line(color='grey') +
  geom_point(size=0.5) + 
  scale_x_discrete(labels = c('pre-IST', 'post-IST')) + 
  scale_y_log10('% T/NK cells', limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100), labels = c('ND', '1', '10','100')) +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
ggsave('./figure_IST/figures/IST1/plots/20230711_mtDNA_dynamics_BM_myeloid.svg', width = 1.1, height = 2)

###

boo = mtDNA.statistics.all.1[, c('mtDNA.mutation', 'TNK_IST1_1_donor', 'TNK_IST1_2_donor')] %>%
  tidyr::pivot_longer(cols = c('TNK_IST1_1_donor', 'TNK_IST1_2_donor'), names_to = 'timepoint')
boo$value[which(boo$value == 0)] = 0.001
ggplot(boo, aes(x=timepoint, y=100*value, group=mtDNA.mutation)) + 
  geom_line(color='grey') +
  geom_point(size=0.5) + 
  scale_x_discrete(labels = c('pre-IST', 'post-IST')) + 
  scale_y_log10('% T/NK cells', limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100), labels = c('ND', '1', '10','100')) +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
ggsave('./figure_IST/figures/IST1/plots/20230711_mtDNA_dynamics_PB_TNK.svg', width = 1.1, height = 2)

boo = mtDNA.statistics.all.1[, c('mtDNA.mutation', 'B cell_IST1_1_donor', 'B cell_IST1_2_donor')] %>%
  tidyr::pivot_longer(cols = c('B cell_IST1_1_donor', 'B cell_IST1_2_donor'), names_to = 'timepoint')
boo$value[which(boo$value == 0)] = 0.001
ggplot(boo, aes(x=timepoint, y=100*value, group=mtDNA.mutation)) + 
  geom_line(color='grey') +
  geom_point(size=0.5) + 
  scale_x_discrete(labels = c('pre-IST', 'post-IST')) + 
  scale_y_log10('% T/NK cells', limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100), labels = c('ND', '1', '10','100')) +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
ggsave('./figure_IST/figures/IST1/plots/20230711_mtDNA_dynamics_PB_Bcell.svg', width = 1.1, height = 2)

boo = mtDNA.statistics.all.1[, c('mtDNA.mutation', 'Myeloid_IST1_1_donor', 'Myeloid_IST1_2_donor')] %>%
  tidyr::pivot_longer(cols = c('Myeloid_IST1_1_donor', 'Myeloid_IST1_2_donor'), names_to = 'timepoint')
boo$value[which(boo$value == 0)] = 0.001
ggplot(boo, aes(x=timepoint, y=100*value, group=mtDNA.mutation)) + 
  geom_line(color='grey') +
  geom_point(size=0.5) + 
  scale_x_discrete(labels = c('pre-IST', 'post-IST')) + 
  scale_y_log10('% T/NK cells', limits = c(0.1, 100), breaks = c(0.1, 1, 10, 100), labels = c('ND', '1', '10','100')) +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
ggsave('./figure_IST/figures/IST1/plots/20230711_mtDNA_dynamics_PB_myeloid.svg', width = 1.1, height = 2)

