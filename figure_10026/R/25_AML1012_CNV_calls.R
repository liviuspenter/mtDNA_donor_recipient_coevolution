# CNV calls for AML1012

library(BuenColors)
library(ComplexHeatmap)
library(data.table)
library(diffloop)
library(stringr)
library(svglite)
library(ArchR)

AML.1012.mito.colors = c('HSCT' = RColorBrewer::brewer.pal(n=5, name='Reds')[5], 
                         'GMP' = RColorBrewer::brewer.pal(n=5, name='Reds')[4], 
                         'Mono' = RColorBrewer::brewer.pal(n=5, name='Reds')[3],
                         'erythroid' = RColorBrewer::brewer.pal(n=5, name='Reds')[2], 
                         'DC' = 'yellow',
                         'CD4' = RColorBrewer::brewer.pal(n=3, name='Blues')[2],
                         'CD8' = RColorBrewer::brewer.pal(n=3, name='Blues')[3],
                         'none' = 'grey')

# cnv calling method provided by Caleb Lareau
source ('./figure_10026/R/cnv_matrix_make_hg38.R')
source ('./R/20200328_mtscatac_seq.R')

source('./figure_10026/R/Tcell.colors.R')

# code adopted from Caleb Lareau
# Quick helper function to grab passed barcodes from cellranger
pull_barcodes_sc_file <- function(sc_file){
  data.frame(fread(sc_file)) %>% filter(cell_id != "None") %>% pull(barcode)
}

# Compute baseline features per bin
AML.1012.mito = loadArchRProject('./data/10026/AML.1012.mito/')

bcs = stringr::str_split_fixed(AML.1012.mito$cellNames[which(AML.1012.mito$Sample == 'AML1012_12' & !AML.1012.mito$manual.clusters %in% c('C4', 'C8'))], pattern = '#', n=2)[,2]
#bcs = bcs[sample(length(bcs), 2000)]
mat = make_basic_mat('./data/10026/AML1012_12/fragments.tsv.gz', bcs)
saveRDS(mat, file = './data/10026/objects/20210612_AML1012_12_cnv.rds')

bcs = stringr::str_split_fixed(AML.1012.mito$cellNames[which(AML.1012.mito$Sample == 'AML1012_34' & !AML.1012.mito$manual.clusters %in% c('C4', 'C8'))], pattern = '#', n=2)[,2]
#bcs = bcs[sample(length(bcs), 2000)]
mat = make_basic_mat('./data/AML1012_34/fragments.tsv.gz', bcs)
saveRDS(mat, file = './data/10026/objects/20210612_AML1012_34_cnv.rds')

base_pbmcs <- readRDS("./data/10026/objects/20210612_Tcell_cnv.rds"); base_pbmcs[is.na(base_pbmcs)] <- 0
cpm_norm <- (t(t(base_pbmcs)/colSums(base_pbmcs)) * 100000)
row_means <- rowMeans(cpm_norm)
row_std <- sqrt(rowVars(cpm_norm))

# Functions to make Z score and log2 change w.r.t. baseline from 10x
makeZscoreMat <- function(mat){
  mat_cpm_norm <- (t(t(mat)/colSums(mat)) * 100000)
  zscore_mat <- (mat_cpm_norm - row_means)/row_std
}

makeLog2Mat <- function(mat){
  mat_cpm_norm <- (t(t(mat)/colSums(mat)) * 100000)
  zscore_mat <- log2(mat_cpm_norm/(row_medians + 1))
  zscore_mat
}

for (s in c('AML1012_1', 'AML1012_2', 'AML1012_4')) {
  if (s %in% c('AML1012_1','AML1012_2')) {
    mat = readRDS('./data/10026/objects/20210612_AML1012_12_cnv.rds')
  } else {
    mat = readRDS('./data/10026/objects/20210612_AML1012_34_cnv.rds')
  }
  
  bcs = c()
  cluster.separators = c()
  
  for (cluster in c('HSCT', 'GMP', 'Mono1', 'Mono2','erythroid', 'DC')) {
    boo = AML.1012.mito$cellNames[which(AML.1012.mito$Sample.hash == s & AML.1012.mito$manual.clusters == cluster)]
    boo = str_split_fixed(boo, pattern = '#', n=2)[,2]
    boo = boo[which(boo %in% colnames(mat))]
    # downsample to 500 cells
    if (length(boo) > 500) {
      boo = boo[sample(length(boo),500)]
    }
    names(boo) = rep(cluster, length(boo))
    bcs = c(bcs, boo)
    cluster.separators = c(cluster.separators, length(bcs))
  }

  mat[is.na(mat)] <- 0
  mat <- mat[,bcs]
  
  # Compute Zscores
  zscore <- makeZscoreMat(mat)
  zscore = zscore[-which(grepl('NA',rownames(zscore))),]
  score_base <- makeZscoreMat(base_pbmcs)
  region_meta <- str_split_fixed(rownames(zscore), "_", 4)
  
  # Cap for visualiation
  zscore[zscore > 3] <- 3
  zscore[zscore < -3] <- -3
  zscore[zscore < 1 & zscore > -1] = 0
  
  keep <- TRUE
  ordering <- factor(region_meta[keep,1], levels = unique(region_meta[keep,1]))
  
  chromosome.separators = c()
  for (i in seq(1,length(unique(ordering)))) {
    chromosome.separators = c(chromosome.separators, sum(table(ordering)[1:i]))
  }
  
  zscore = t(zscore)
  zscore = zscore[bcs,]
  
  side.colors = c(AML.1012.mito.colors[names(bcs)])
  
  png(paste0('./figure_10026/figures/CNV/20210612_',s,'_cnv.png'), width = 4*1200, height = 1*1200, res = 600)
  heatmap.2(zscore, dendrogram = 'none', Rowv = 'True', Colv = 'None', key = F, 
            col = as.character(jdb_palette("solar_basic",type="continuous")), 
            breaks = seq(-3,3,0.006), 
            trace = 'none', scale = 'none', labRow = F, labCol = F,
            colsep = chromosome.separators, 
            rowsep = cluster.separators,
            RowSideColors = side.colors,
            margins = c(0,0), lhei = c(0.01,1), lwid = c(0.1,10),
            sepcolor = 'white', sepwidth = c(0.05,0.05))
  dev.off()
  
  message(paste0(s, ' ', length(bcs), ' cells'))
}

### focus on del(5q)(q11.2q34) and amp(22q)
### q11.2 starts at 51,400,001
### q34 starts at 160,500,001

### 22q starts at 15,000,001

zscore.all = matrix()
for (s in c('AML1012_1', 'AML1012_2', 'AML1012_4')) {
  message(s)
  if (s %in% c('AML1012_1','AML1012_2')) {
    mat = readRDS('./data/10026/objects/20210612_AML1012_12_cnv.rds')
  } else {
    mat = readRDS('./data/10026/objects/20210612_AML1012_34_cnv.rds')
  }
  
  bcs = c()
  cluster.separators = c()
  
  for (cluster in c('HSCT', 'GMP', 'Mono','erythroid', 'DC', 'CD4', 'CD8')) {
    boo = AML.1012.mito$cellNames[which(AML.1012.mito$Sample.hash == s & AML.1012.mito$manual.clusters == cluster)]
    boo = str_split_fixed(boo, pattern = '#', n=2)[,2]
    boo = boo[which(boo %in% colnames(mat))]
    names(boo) = rep(cluster, length(boo))
    bcs = c(bcs, boo)
    cluster.separators = c(cluster.separators, length(bcs))
  }
  
  mat[is.na(mat)] <- 0
  mat <- mat[,bcs]
  
  # Compute Zscores
  zscore <- makeZscoreMat(mat)
  zscore = zscore[-which(grepl('NA',rownames(zscore))),]
  score_base <- makeZscoreMat(base_pbmcs)
  region_meta <- str_split_fixed(rownames(zscore), "_", 4)
  
  # Cap for visualiation
  #zscore[zscore > 3] <- 3
  #zscore[zscore < -3] <- -3
  #zscore[zscore < 1 & zscore > -1] = 0
  
  if (s %in% c('AML1012_1', 'AML1012_2')) {
    colnames(zscore) = paste0('AML1012_12#', colnames(zscore))
  } else {
    colnames(zscore) = paste0('AML1012_34#', colnames(zscore))
  }
  
  
  if (s == 'AML1012_1') {
    zscore.all = zscore
  } else {
    zscore.all = cbind(zscore.all, zscore)
  }
}
boo = data.frame(bc = AML.1012.mito$cellNames, cluster = AML.1012.mito$manual.clusters, sample = AML.1012.mito$Sample.hash)
rownames(boo) = boo$bc
#boo = boo[which(boo$bc %in% rownames(del.5q)),]

regions.5q = rownames(zscore)[grepl('chr5', rownames(zscore))]
regions.5q = regions.5q[seq(which(regions.5q == 'chr5_52000000_62000000_27'), which(regions.5q == 'chr5_150000000_160000000_76'))]

regions.22q = rownames(zscore)[grepl('chr22', rownames(zscore))]
regions.22q = regions.22q[seq(which(regions.22q == 'chr22_16000000_26000000_9'), which(regions.22q == 'chr22_40000000_50000000_21'))]

CNV.df = data.frame(del.5q = as.numeric(colMeans(zscore.all[regions.5q,])),
                    amp.22q = as.numeric(colMeans(zscore.all[regions.22q,])))
rownames(CNV.df) = colnames(zscore.all)
CNV.df$celltype = boo[rownames(CNV.df), 'cluster']
CNV.df$sample = boo[rownames(CNV.df), 'sample']

ggplot(data=CNV.df, aes(x=del.5q, y=amp.22q)) + 
  geom_density_2d(bins=100, color='black') + 
  geom_point(aes(color=celltype), size=0.5) +
  scale_x_continuous('del(5q)',limits = c(-2,2)) +
  scale_y_continuous('amp(22q)',limits = c(-2,3)) +
  scale_color_manual(values = AML.1012.mito.colors) + 
  geom_vline(xintercept = -0.5) +
  geom_hline(yintercept = 1) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figure_10026/figures/CNV/20220302_AML1012_CNV_plot.svg', width = 2, height = 2)

ggplot(CNV.df, aes(x=del.5q)) + geom_histogram(binwidth = 0.1) + geom_vline(xintercept = -0.5) + facet_wrap(~celltype)
ggplot(CNV.df, aes(x=amp.22q)) + geom_histogram(binwidth = 0.1) + geom_vline(xintercept = 1) + facet_wrap(~celltype)

CNV.df$del.5q.detected = ifelse(CNV.df$del.5q < -0.5, 1, 0)
CNV.df$amp.22q.detected = ifelse(CNV.df$amp.22q > 1, 1, 0)
CNV.df$normal = ifelse(CNV.df$del.5q.detected == 0 & CNV.df$amp.22q.detected == 0, 1, 0)

df = getEmbedding(AML.1012.mito)
colnames(df) = c('UMAP1', 'UMAP2')
df$CNV = 'none'
df[rownames(CNV.df),'normal'] = CNV.df$normal
df[AML.1012.mito$cellNames,'celltype'] = AML.1012.mito$manual.clusters
df[AML.1012.mito$cellNames,'Sample'] = AML.1012.mito$Sample.hash

ggplot() + 
  geom_point(data=df[which(df$normal == 0),], aes(x=UMAP1, y=UMAP2), color='lightgrey', size=0.5) +
  geom_point(data=df[which(df$normal == 1 & df$celltype %in% c('CD4', 'CD8')),], aes(x=UMAP1, y=UMAP2), color='lightblue', size=1) +
  geom_point(data=df[which(df$normal == 1 & (!df$celltype %in% c('CD4', 'CD8'))),], aes(x=UMAP1, y=UMAP2), color='black', size=1) +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
ggsave('./figure_10026/figures/umaps/20220301_AML1012_CNV_cells.png', width = 2.5, height = 2.5)

length(which(df$celltype %in% c('CD4', 'CD8')))
length(which(df$celltype %in% c('HSCT', 'GMP', 'Mono1', 'Mono2','erythroid', 'DC') & df$normal == 0))
length(which(df$celltype %in% c('HSCT', 'GMP', 'Mono1', 'Mono2','erythroid', 'DC') & df$normal == 1))

boo = CNV.df %>% filter(celltype %in% c('HSCT', 'GMP', 'Mono', 'erythroid')) %>% 
  group_by(sample, celltype) %>% 
  summarize(normal = length(which(normal == 1)) / length(which(normal == 0)))

ggplot(boo, aes(x=sample, y=100*normal, color=celltype)) + geom_line(aes(group=celltype)) + geom_point(size=0.5) +
  scale_x_discrete(labels = c('Screening', 'EOLN', 'EOT')) + 
  scale_y_continuous('% cells without CNV') +
  scale_color_manual(values = AML.1012.mito.colors) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank())
ggsave('./figure_10026/figures/CNV/20220302_AML1012_CNV_longitudinal.svg', width = 1, height = 2)
