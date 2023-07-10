# process mtDNA mutations to create list of maternal mtDNA variants for each patient
# create ArchR object with mtDNA coverage >10x

library(ArchR)
library(ggplot2)
library(ggrepel)
library(gplots)
library(parallel)
library(Signac)
library(Seurat)
library(BuenColors)
library(dplyr)

# combine mtDNA mutations
combined.mutation.frequencies.12 = readRDS(file = './data/IST/mtDNA/20220110_IST1_2_combined_mutation_frequencies.rds')
combined.mutation.frequencies.12$mutation = rownames(combined.mutation.frequencies.12)
combined.mutation.frequencies.3 = readRDS(file = './data/IST/mtDNA/20220110_IST3_combined_mutation_frequencies.rds')
combined.mutation.frequencies.3$mutation = rownames(combined.mutation.frequencies.3)
combined.mutation.frequencies.4 = readRDS(file = './data/IST/mtDNA/20220110_IST4_combined_mutation_frequencies.rds')
combined.mutation.frequencies.4$mutation = rownames(combined.mutation.frequencies.4)
combined.mutation.frequencies.5 = readRDS(file = './data/IST/mtDNA/20220110_IST5_combined_mutation_frequencies.rds')
combined.mutation.frequencies.5$mutation = rownames(combined.mutation.frequencies.5)

combined.mutation.frequencies = merge(combined.mutation.frequencies.12, combined.mutation.frequencies.3, by='mutation', all=T)
combined.mutation.frequencies = merge(combined.mutation.frequencies, combined.mutation.frequencies.4, by='mutation', all=T)
combined.mutation.frequencies = merge(combined.mutation.frequencies, combined.mutation.frequencies.5, by='mutation', all=T)

rownames(combined.mutation.frequencies) = combined.mutation.frequencies$mutation
combined.mutation.frequencies = combined.mutation.frequencies[,-1]
saveRDS(combined.mutation.frequencies, file = './data/IST/mtDNA/20220110_IST_all_combined_mutation_frequencies.rds')

# extract donor and recipient variants
combined.mutation.frequencies.12 = readRDS(file = './data/IST/mtDNA/20220110_IST1_2_combined_mutation_frequencies.rds')
combined.mutation.frequencies.3 = readRDS(file = './data/IST/mtDNA/20220110_IST3_combined_mutation_frequencies.rds')
combined.mutation.frequencies.4 = readRDS(file = './data/IST/mtDNA/20220110_IST4_combined_mutation_frequencies.rds')
combined.mutation.frequencies.5 = readRDS(file = './data/IST/mtDNA/20220110_IST5_combined_mutation_frequencies.rds')
germline.variants = data.frame()
# IST1 (samples IST1 + IST2)
ComplexHeatmap::Heatmap(combined.mutation.frequencies.12[which(rowMeans(combined.mutation.frequencies.12) > 0.01),
                                                         sample(ncol(combined.mutation.frequencies.12), size=1000)], show_column_names = F)

donor.variants = c('16224T>C', '1189T>C', '5237G>A', '1811A>G', '14798T>C', '10154A>G', '12738T>G', '11299T>C', '9896A>G', '5913G>A',
                   '10398A>G', '16311T>C', '10550A>G', '6845C>T', '15301G>A', '146T>C', '16093T>C')
recipient.variants = c('7705T>C', '10007T>C', '16324T>C', '14323G>A', '15829A>G', '7808C>T', '5165C>T', '16290C>T', '16234C>T', 
                       '16221C>T', '16189T>C', '16183A>C')

germline.variants = rbind(germline.variants, data.frame(
  variant = c(donor.variants, recipient.variants),
  individual = c(rep('donor', length(donor.variants)), rep('recipient', length(recipient.variants))),
  sample = rep('IST1', length(c(donor.variants, recipient.variants)))
))

# IST2 (samples IST3)
ComplexHeatmap::Heatmap(combined.mutation.frequencies.3[which(rowMeans(combined.mutation.frequencies.3) > 0.01),
                                                         sample(ncol(combined.mutation.frequencies.3), size=1000)], show_column_names = F)

donor.variants = c('13953T>C', '11253T>C', '16362T>C', '3915G>A', '9380G>A', '4727A>G', '239T>C', '16482A>G')
recipient.variants = c('4216T>C', '16519T>C', '16297T>C', '10463T>C', '73A>G', '15452C>A', '16126T>C', '13887A>G', '15607A>G', 
                       '13368G>A', '4917A>G', '5147G>A', '8697G>A', '11812A>G', '12441T>C', '2706A>G', '7028C>T', '14766C>T', 
                       '1888G>A', '14233A>G', '15222A>G', '709G>A', '16296C>T', '7984G>A', '11251A>G', '7354T>C',
                       '11719G>A', '16304T>C', '16294C>T', '930G>A', '14905G>A', '15928G>A', '9224T>C', '2885T>C',
                       '16147C>T', '152T>C')

germline.variants = rbind(germline.variants, data.frame(
  variant = c(donor.variants, recipient.variants),
  individual = c(rep('donor', length(donor.variants)), rep('recipient', length(recipient.variants))),
  sample = rep('IST3', length(c(donor.variants, recipient.variants)))
))

# IST3 (samples IST4)
ComplexHeatmap::Heatmap(combined.mutation.frequencies.4[which(rowMeans(combined.mutation.frequencies.4) > 0.01),
                                                        sample(ncol(combined.mutation.frequencies.4), size=1000)], show_column_names = F)

donor.variants = c('200A>G','189A>G', '16292C>T', '16311T>C', '10070C>T', '1822T>C', '3396T>C', '9950T>C', '3505A>G', '14769A>G', '16209T>C', '16218C>T',
                   '5601C>T', '8527A>G', '11440G>A', '15514T>C', '7819C>A', '8932C>T', '16519T>C', '4218T>C', '150C>T')
recipient.variants = c('16093T>C', '146T>C', '152T>C', '12693A>G', '10115T>C', '2416T>C', '14566A>G', '13803A>G', '3594C>T', '15784T>C', '7771A>G',
                       '11944T>C', '11914G>A', '4104A>G', '13590G>A', '7175T>C', '2789C>T', '1018G>A', '769G>A', '5319A>G', '13650C>T', '7521G>A',
                       '9221A>G', '5581A>G', '8206G>A', '7256C>T', '7274C>T', '16390G>A', '16278C>T', '16294C>T', '16309A>G', '16189T>C', '195T>C')

germline.variants = rbind(germline.variants, data.frame(
  variant = c(donor.variants, recipient.variants),
  individual = c(rep('donor', length(donor.variants)), rep('recipient', length(recipient.variants))),
  sample = rep('IST4', length(c(donor.variants, recipient.variants)))
))

# IST4 (samples IST5)
ComplexHeatmap::Heatmap(combined.mutation.frequencies.5[which(rowMeans(combined.mutation.frequencies.5) > 0.01),], show_column_names = F)

donor.variants = c('16304T>C', '16294C>T', '16296C>T', '15607A>G', '11251A>G', '2706A>G', '13368G>A', '930G>A', '1888G>A', '15928G>A', '709G>A',
                   '14905G>A', '7028C>T', '14766C>T', '13928G>C', '4216T>C', '4917A>G', '11719G>A', '11812A>G', '15452C>A', '8697G>A', '5147G>A',
                   '14233A>G', '73A>G')
recipient.variants = c('4793A>G', '10304T>C', '3984C>T', '11016G>A', '11971C>T')

germline.variants = rbind(germline.variants, data.frame(
  variant = c(donor.variants, recipient.variants),
  individual = c(rep('donor', length(donor.variants)), rep('recipient', length(recipient.variants))),
  sample = rep('IST5', length(c(donor.variants, recipient.variants)))
))

write.csv2(germline.variants, file = './data/IST/objects/20220117_IST_germline_variants.csv', quote = F, row.names = F)
germline.variants = data.table::fread(file = './data/IST/objects/20220117_IST_germline_variants.csv')

# subset ArchR archive to include only cells with mtDNA coverage >10x
IST.asap = loadArchRProject('./data/IST/IST.asap/') 

IST.asap.mito = subsetArchRProject(ArchRProj = IST.asap, 
                                   cells = IST.asap$cellNames[which(IST.asap$cellNames %in% colnames(combined.mutation.frequencies))],
                                   outputDirectory = './data/IST/IST.asap.mito/', threads = 12, force = T)
IST.asap.mito = addClusters(input = IST.asap.mito, reducedDims = 'IterativeLSI', method = 'Seurat', name = 'Clusters', resolution = 0.3, force = T)
IST.asap.mito = addUMAP(ArchRProj = IST.asap.mito, reducedDims = 'IterativeLSI', name = 'UMAP', nNeighbors = 30, minDist = 0.5, metric = 'cosine', force = T)
IST.asap.mito = addImputeWeights(IST.asap.mito)

IST.asap.mito <- addGroupCoverages(IST.asap.mito, groupBy = "manual.cluster", threads = 12)
IST.asap.mito = addReproduciblePeakSet(IST.asap.mito, groupBy = 'manual.cluster', threads = 12)
IST.asap.mito = addPeakMatrix(IST.asap.mito, threads = 12)
IST.asap.mito <- addMotifAnnotations(IST.asap.mito, motifSet = "cisbp", name = "Motif")
IST.asap.mito <- addBgdPeaks(IST.asap.mito)
IST.asap.mito <- addDeviationsMatrix(IST.asap.mito, peakAnnotation = "Motif", threads = 12, force = T)

saveArchRProject(IST.asap.mito)
IST.asap.mito = loadArchRProject('./data/IST/IST.asap.mito/')

IST.asap.mito$patient = stringr::str_split_fixed(IST.asap.mito$Sample, pattern = '_', n=2)[,1]
IST.asap.mito$patient[which(IST.asap.mito$patient == 'IST2')] = 'IST1'

IST.asap.mito$individual = 'none'

df.all = data.frame()
for (patient in unique(germline.variants$sample)) {
  cells = IST.asap.mito$cellNames[which(IST.asap.mito$patient == patient)]
  
  df = data.frame(donor = colMeans(combined.mutation.frequencies[germline.variants$variant[which(germline.variants$sample == patient & 
                                                                                                   germline.variants$individual == 'donor')],cells]),
                  recipient = colMeans(combined.mutation.frequencies[germline.variants$variant[which(germline.variants$sample == patient & 
                                                                                                   germline.variants$individual == 'recipient')],cells]))
  
  df$individual = 'none'
  df$individual[which(df$donor > 0.8 & df$recipient < 0.2)] = 'donor'
  df$individual[which(df$donor < 0.2 & df$recipient > 0.8)] = 'recipient'
  rownames(df) = cells
  df.all = rbind(df.all, df)
    
  
  p=ggplot(df, aes(x=100*recipient, y=100*donor, color=individual)) + geom_point(size=0.5) + 
    scale_x_continuous('% recipient') + 
    scale_y_continuous('% donor') + 
    scale_color_manual(values = c('none' = 'grey', 'donor' = 'purple', 'recipient' = 'orange')) + 
    geom_hline(yintercept = c(20,80)) + 
    geom_vline(xintercept = c(20,80)) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title = element_text('Arial', size=10, color='black'))
  ggsave(paste0('./figure_IST/figures/deconvolution/20220117_',patient,'_mtDNA_variants.svg'), width = 1.5, height = 1.5, plot = p)
}
IST.asap.mito$individual = df.all[IST.asap.mito$cellNames, 'individual']

saveArchRProject(IST.asap.mito)