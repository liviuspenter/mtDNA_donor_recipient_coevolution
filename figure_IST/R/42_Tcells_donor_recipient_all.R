# donor-recipient CD4/CD8 T cells

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

donor.cells = rownames(IST.asap.mito)[which(IST.asap.mito$manual.cluster == 'T' & IST.asap.mito$individual == 'donor')]
recipient.cells = rownames(IST.asap.mito)[which(IST.asap.mito$manual.cluster == 'T' & IST.asap.mito$individual == 'recipient')]

TSB.so = readRDS('./data/IST/objects/20220118_IST_mito_TSB.so.rds')

TSB.so$group = 'none'
TSB.so$group[donor.cells] = 'donor'
TSB.so$group[recipient.cells] = 'recipient'
TSB.so$group = factor(TSB.so$group, levels = c('donor', 'recipient', 'none'))

df = GetAssayData(TSB.so, slot = 'scale.data')[c('CD4', 'CD8'), c(donor.cells, recipient.cells)] %>%
  t %>% as.data.frame()

df$group = TSB.so$group[rownames(df)]

ggplot(df, aes(x=CD4, y=CD8)) + geom_point(aes(color=group), size=0.5) +
  scale_color_manual(values = c('recipient' = 'purple', 'donor' = 'orange')) + 
  geom_hline(yintercept = 0.5) + 
  geom_vline(xintercept = 0.5) + 
  theme_classic() + 
  theme(legend.position = 'none',
        plot.title = element_blank(),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_IST/figures/plots/20230424_CD4_CD8_donor_recipient.svg', width = 1.5, height = 1.5)

# statistics
df$phenotype = ifelse(df$CD4 > 0.5 & df$CD8 < 0.5, 'CD4', 'other')
df$phenotype[which(df$CD4 < 0.5 & df$CD8 > 0.5)] = 'CD8'

boo = df %>% group_by(group) %>% summarize(CD4 = length(which(phenotype == 'CD4')),
                                     CD8 = length(which(phenotype == 'CD8')))

boo$CD4 / (boo$CD4 + boo$CD8)
fisher.test(as.matrix(boo[,c(2,3)]))
