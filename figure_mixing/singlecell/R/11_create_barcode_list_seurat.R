library(ArchR)
library(ggplot2)
library(Seurat)

CLL.mix.so = readRDS('./data/mixing/singlecell/CLL_mix_so.rds')

# barcodes of both libraries
file1 = colnames(subset(CLL.mix.so, orig.ident == 'CLL1'))
file2 = colnames(subset(CLL.mix.so, orig.ident == 'CLL3'))
file1 = stringr::str_split_fixed(file1, pattern = '_', n=2)[,2]
file2 = stringr::str_split_fixed(file2, pattern = '_', n=2)[,2]

# remove barcodes that overlap between both libraries
file1 = setdiff(file1, file2)
file2 = setdiff(file2, file1)

for (downsample in c(1,5,10,50,100,500,1000)) {
  file1.downsample = as.data.frame(file1[sample(length(file1), downsample)])
  file2.downsample = as.data.frame(file2[sample(length(file2), 10000)])
  file1.downsample$V2 = 'mixing'
  file2.downsample$V2 = 'mixing'
  write.table(file1.downsample, file = paste0('./data/mixing/singlecell/barcodes_scRNAseq/CLL_relapse1_1_barcodes.',downsample), col.names = F, row.names = F, quote = F)
  write.table(file2.downsample, file = paste0('./data/mixing/singlecell/barcodes_scRNAseq/CLL_relapse3_1_barcodes.',downsample), col.names = F, row.names = F, quote = F)
}