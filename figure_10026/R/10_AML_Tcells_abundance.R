# T cell subset dynamics

library(ArchR)
library(gplots)
library(Seurat)

source('./figure_10026/R/Tcell.colors.R')

AML.Tcell = loadArchRProject('./data/10026/AML.Tcell/')
TSB.Tcell = readRDS('./data/10026/objects/20210603_TSB_Tcell.rds')

cluster.mat = data.frame(cluster = as.character(), sample = as.character(), cells = as.numeric())
for (cluster in unique(AML.Tcell$manual_clusters)) {
  for (sample in unique(AML.Tcell$Sample.hash)) {
    cluster.mat = rbind(cluster.mat, data.frame(cluster = cluster, 
                                                sample = sample, 
                                                cells = length(which(AML.Tcell$manual_clusters == cluster & 
                                                                       AML.Tcell$Sample.hash == sample))))
  }
}
cluster.mat$cells.freq = 0
for (i in seq(1,nrow(cluster.mat))) {
  cluster.mat$cells.freq[i] = cluster.mat$cells[i] / sum(cluster.mat$cells[which(cluster.mat$sample == cluster.mat$sample[i])])
}
cluster.mat$cluster = factor(cluster.mat$cluster, levels = names(Tcell.colors))
p=ggplot(cluster.mat[which(cluster.mat$sample != 'none'),], aes(x=sample, y=cells.freq, fill=cluster)) + 
  geom_col() + 
  scale_fill_manual(values = Tcell.colors) + 
  scale_x_discrete(labels = c('Screening', 'EOLN', 'C1', 'C4', 'C10', 'Screening', 'EOLN', 'EOT', 'Screening', 'EOLN', 'C1', 'Relapse')) + 
  scale_y_continuous('%cells',labels = scales::percent) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank(),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave(filename = './figure_10026/figures/plots/20210604_Tcell_cluster.svg', width = 2.5, height = 3, plot = p)