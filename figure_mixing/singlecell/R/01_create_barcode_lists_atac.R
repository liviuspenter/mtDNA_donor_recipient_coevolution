# generate barcodes for each mixing step for each sample
# cell barcodes that overlap between samples are excluded

library(ArchR)
library(dplyr)
library(ggplot2)
library(parallel)
library(Seurat)
library(Signac)

source('./R/variantplot.LP.R')

CLL.mix = loadArchRProject('./data/mixing/singlecell/CLL.mix/')

# read barcode files of both libraries
file1 = data.table::fread('./data/mixing/singlecell/barcodes_mtscATACseq/CLL_relapse1_1_barcodes.tsv.gz', header = F)
file2 = data.table::fread('./data/mixing/singlecell/barcodes_mtscATACseq/CLL_relapse3_1_barcodes.tsv.gz', header = F)

# remove barcodes that overlap between both libraries
file1 = file1 %>% filter(!V1 %in% file2$V1)
file2 = file2 %>% filter(!V1 %in% file1$V1)

# remove barcodes that are not in the ArchR object
file1 = file1 %>% filter(paste0('CLL_relapse1_1#',V1) %in% rownames(CLL.mix))
file2 = file2 %>% filter(paste0('CLL_relapse3_1#',V1) %in% rownames(CLL.mix))

# write barcodes
for (downsample in c(1,5,10,50,100,500,1000)) {
  file1.downsample = file1[sample(nrow(file1), downsample),]
  file2.downsample = file2#[sample(nrow(file2), 8000),]
  file1.downsample$V2 = 'mixing'
  file2.downsample$V2 = 'mixing'
  write.table(file1.downsample, file = paste0('./data/mixing/singlecell/barcodes_mtscATACseq/CLL_relapse1_1_barcodes.',downsample), col.names = F, row.names = F, quote = F)
  write.table(file2.downsample, file = paste0('./data/mixing/singlecell/barcodes_mtscATACseq/CLL_relapse3_1_barcodes.',downsample), col.names = F, row.names = F, quote = F)
}
