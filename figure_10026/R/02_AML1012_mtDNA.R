# process mgatk output of AML1012

library(ggplot2)
library(ggrepel)
library(gplots)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(BuenColors)
library(dplyr)

samples = c('AML1012_12','AML1012_34')
atac.paths = c('./data/10026/AML1012_12/',
               './data/10026/AML1012_34/')
mgatk.paths = c('./data/10026/AML1012_12.mgatk/',
               './data/10026/AML1012_34.mgatk/')

# load gene annotations from Ensembl
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

min_coverage = 20

for (i in seq(1,length(samples))) {
  message(paste0('Extracing barcodes for sample ',samples[i]))
  
  counts <- Read10X_h5(filename = paste0(atac.paths[i], 'filtered_peak_bc_matrix.h5'))
  metadata <- read.csv(file = paste0(atac.paths[i], 'singlecell.csv'), header = TRUE, row.names = 1)
  
  # create object
  crc_assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), annotation = annotations, min.cells = 10, genome = "hg38",
    fragments = paste0(atac.paths[i], 'fragments.tsv.gz'))
  crc <- CreateSeuratObject(counts = crc_assay, assay = 'peaks', meta.data = metadata)
  
  # Augment QC metrics that were computed by cellranger-atac
  crc$pct_reads_in_peaks <- crc$peak_region_fragments / crc$passed_filters * 100
  crc$pct_reads_in_DNase <- crc$DNase_sensitive_region_fragments / crc$passed_filters * 100
  crc$blacklist_ratio <- crc$blacklist_region_fragments / crc$peak_region_fragments
  
  # compute TSS enrichment score and nucleosome banding pattern
  crc <- TSSEnrichment(crc)
  crc <- NucleosomeSignal(crc)
  
  # remove low-quality cells
  crc <- subset(
    x = crc,
    subset = nCount_peaks > 1000 &
      nCount_peaks < 50000 &
      pct_reads_in_DNase > 40 &
      blacklist_ratio < 0.05 &
      #TSS.enrichment > 3 &
      nucleosome_signal < 4
  )
  
  # load mgatk output and process 
  mito.data <- ReadMGATK(dir = mgatk.paths[i])
  mito <- CreateAssayObject(counts = mito.data$counts)
  mito <- subset(mito, cells = colnames(crc))
  crc[["mito"]] <- mito
  crc <- AddMetaData(crc, metadata = mito.data$depth, col.name = "mtDNA_depth")
  
  p=qplot(crc$mtDNA_depth) + scale_x_log10(breaks = c(1,10,100)) +
    geom_vline(xintercept = 10, color = "firebrick") +
    pretty_plot() + L_border() + labs(x = "mean mtDNA depth/bp", y = "count")
  ggsave(plot = , p,filename = paste0('./figure_10026/figures/qc/20210429_',samples[i],'_qc_plot1.png'), 
         width = 4, height = 4, dpi = 350)
  
  crc <- subset(crc, mtDNA_depth >= min_coverage)
  
  # save Seurat object and barcodes
  saveRDS(object = crc, file = paste0('./data/10026/objects/20210429_',samples[i],'_seuratobject.rds'))
  write.table(data.frame(colnames(crc)), 
              file = paste0('./data/10026/barcodes/',samples[i],'_barcodes.csv'), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# combine mitochondrial mutations from multiple samples to identify variants based on all samples together
for (i in seq(1,length(samples))) {
  message(paste0('Reading sample ',samples[i]))
  
  data.mgatk = ReadMGATK(dir = paste0(mgatk.paths[i]))
  barcodes = read.table(file = paste0('./data/barcodes/',samples[i],'_barcodes.csv'))
  quality.cells = which(colnames(data.mgatk$counts) %in% barcodes$V1)
  data.mgatk$counts = data.mgatk$counts[,quality.cells]
  data.mgatk$depth = subset(data.mgatk$depth, rownames(data.mgatk$depth) %in% colnames(data.mgatk$counts))
  colnames(data.mgatk$counts) = paste0(samples[i],'#',colnames(data.mgatk$counts))
  rownames(data.mgatk$depth) = paste0(samples[i],'#',rownames(data.mgatk$depth))
  
  if (i == 1) {
    combined.data = data.mgatk
  } else {
    combined.data$counts = cbind(combined.data$counts, data.mgatk$counts)
    combined.data$depth = rbind(combined.data$depth, data.mgatk$depth)
  }
}
combind.mutations = IdentifyVariants(object = combined.data$counts, refallele = combined.data$refallele)
p=VariantPlot(variants = combind.mutations)
ggsave(filename = paste0('./figure_10026/figures/qc/20210429_AML1012_strand_concordance_vmr.png'), plot = p, 
       width = 4, height = 4)

# Establish a filtered data frame of variants based on this processing
high.conf <- subset(combind.mutations, subset = n_cells_conf_detected >= 5 & strand_correlation >= 0.65 & vmr > 0.01)
write.table(high.conf, file = './data/10026/mtDNA/20210429_AML1012_high_confidence_variants.csv', sep = '\t', row.names = F, quote = F)

# calculate frequencies for all cells
combined.mutation.frequencies = data.frame(mutation = as.character())
for (i in seq(1,length(samples))) {
  crc = readRDS(file = paste0('./data/10026/objects/20210429_',samples[i],'_seuratobject.rds'))
  
  crc <- AlleleFreq(object = crc, variants = high.conf$variant, assay = "mito")
  boo=as.data.frame(GetAssayData(crc, assay = 'alleles'))
  colnames(boo) = paste0(samples[i],'#',colnames(boo))
  boo$mutation = rownames(boo)
  
  combined.mutation.frequencies=full_join(combined.mutation.frequencies, boo, by = 'mutation')
  saveRDS(object = crc, file = paste0('./data/10026/objects/20210429_',samples[i],'_seuratobject.rds'))
}
rownames(combined.mutation.frequencies) = combined.mutation.frequencies$mutation
combined.mutation.frequencies = combined.mutation.frequencies[,-1]
saveRDS(combined.mutation.frequencies, file = './data/10026/mtDNA/20210429_AML1012_combined_mutation_frequencies.rds')
