setwd('/Users/shaka87/dfci/mito_atac/')

library(BuenColors)
library(ComplexHeatmap)
library(data.table)
library(stringr)
library(svglite)
library(ArchR)

# cnv calling method provided by Caleb Lareau
source ('./data/cnv_cl/00_cnv_matrix_make_hg38.R')
source('./data/cnv_cl/20210529_cnv_normal.R')
source ('/Users/shaka87/dfci/scripts/20200328_mtscatac_seq.R')

setwd('/Users/shaka87/dfci/asap_seq/')
source('./analysis/Tcell.colors.R')

# code adopted from Caleb Lareau
# Quick helper function to grab passed barcodes from cellranger
pull_barcodes_sc_file <- function(sc_file){
  data.frame(fread(sc_file)) %>% filter(cell_id != "None") %>% pull(barcode)
}

# Compute baseline features per bin
AML.1026.mito = loadArchRProject('/Users/shaka87/ArchRProjects/AML.asap/AML.1026.mito/')

bcs = stringr::str_split_fixed(AML.1026.mito$cellNames[which(AML.1026.mito$Sample == 'AML1026_12' & !AML.1026.mito$manual.clusters %in% c('T cell'))], pattern = '#', n=2)[,2]
bcs = bcs[sample(length(bcs), 2000)]
mat = make_basic_mat('./data/AML1026_12/fragments.tsv.gz', bcs)
saveRDS(mat, file = './data/objects/20210612_AML1026_12_cnv.rds')

bcs = stringr::str_split_fixed(AML.1026.mito$cellNames[which(AML.1026.mito$Sample == 'AML1026_34' & !AML.1026.mito$manual.clusters %in% c('T cell'))], pattern = '#', n=2)[,2]
bcs = bcs[sample(length(bcs), 2000)]
mat = make_basic_mat('./data/AML1026_34/fragments.tsv.gz', bcs)
saveRDS(mat, file = './data/objects/20210612_AML1026_34_cnv.rds')

base_pbmcs <- readRDS("./data/objects/20210612_Tcell_cnv.rds"); base_pbmcs[is.na(base_pbmcs)] <- 0
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

for (s in c('AML1026_1', 'AML1026_2', 'AML1026_3','AML1026_4')) {
  if (s %in% c('AML1026_1', 'AML1026_2')) {
    mat = readRDS('./data/objects/20210612_AML1026_12_cnv.rds')
  } else {
    mat = readRDS('./data/objects/20210612_AML1026_34_cnv.rds')
  }
  
  bcs = c()
  for (cluster in c('HSCT', 'GMP', 'Mono', 'erythroid')) {
    boo = AML.1026.mito$cellNames[which(AML.1026.mito$Sample.hash == s & AML.1026.mito$manual.clusters == cluster)]
    boo = str_split_fixed(boo, pattern = '#', n=2)[,2]
    boo = boo[which(boo %in% colnames(mat))]
    # downsample to 500 cells
    if (length(boo) > 500) {
      boo = boo[sample(length(boo),500)]
    }
    names(boo) = rep(cluster, length(boo))
    bcs = c(bcs, boo)
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
  
  side.colors = c(AML.1026.mito.colors[names(bcs)])
  
  png(paste0('./figures/cnv/20210612_',s,'_cnv.png'), width = 4*1200, height = 1*1200, res = 600)
  heatmap.2(zscore, dendrogram = 'none', Rowv = 'True', Colv = 'None', key = F, 
            col = as.character(jdb_palette("solar_basic",type="continuous")), 
            breaks = seq(-3,3,0.006), 
            trace = 'none', scale = 'none', labRow = F, labCol = F,
            colsep = chromosome.separators, 
            #rowsep = cluster.separators,
            RowSideColors = side.colors,
            margins = c(0,0), lhei = c(0.01,1), lwid = c(0.1,10),
            sepcolor = 'white', sepwidth = c(0.05,0.05))
  dev.off()
}
