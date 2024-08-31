library(ArchR)
library(dplyr)
library(ggplot2)
library(parallel)
library(Seurat)
library(Signac)

source("./R/variantplot.LP.R")

AML.mix <- loadArchRProject("./data/mixing/AML.mix/")

# read barcode files of both libraries
file1 <- data.table::fread("./data/10026/Pool148_13/filtered_peak_bc_matrix/barcodes.tsv", header = F)
file2 <- data.table::fread("./data/10026/Pool148_17/filtered_peak_bc_matrix/barcodes.tsv", header = F)

# remove barcodes that overlap between both libraries
file1 <- file1 %>% filter(!V1 %in% file2$V1)
file2 <- file2 %>% filter(!V1 %in% file1$V1)

# remove barcodes that are not in the ArchR object
file1 <- file1 %>% filter(paste0("AML1011_A#", V1) %in% AML.mix$cellNames)
file2 <- file2 %>% filter(paste0("AML1012_A#", V1) %in% AML.mix$cellNames)

for (downsample in c(1, 5, 10, 50, 100, 500, 1000)) {
  file1.downsample <- file1[sample(nrow(file1), downsample), ]
  file2.downsample <- file2 # [sample(nrow(file2), 8000),]
  file1.downsample$V2 <- "mixing"
  file2.downsample$V2 <- "mixing"
  write.table(file1.downsample, file = paste0("./data/mixing/singlecell/artificial_mixing_AML/barcodes_mtscATACseq/AML1011_barcodes.", downsample), col.names = F, row.names = F, quote = F)
  write.table(file2.downsample, file = paste0("./data/mixing/singlecell/artificial_mixing_AML/barcodes_mtscATACseq/AML1012_barcodes.", downsample), col.names = F, row.names = F, quote = F)
}
