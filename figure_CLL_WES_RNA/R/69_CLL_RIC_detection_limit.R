library(cowplot)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(scales)
library(tidyr)

meta = readxl::read_excel('./data/CLL_RIC/20231229_CLL_RIC_mtDNA_meta.xlsx') %>%
  as.data.frame()
meta2 = readxl::read_excel('./data/CLL_RIC/20231229_CLL_RIC_mtDNA_meta.xlsx', sheet = 2) %>%
  as.data.frame()
rownames(meta2) = meta2$Patient
meta$Day2 = 0
meta[which(meta$Timepoint == 2), 'Day2'] = meta2[as.character(meta$Patient[which(meta$Timepoint == 2)]), 'Days_between_samples']


### read all data
mgatk.data.RNA = 0 
for (library in list.files('./data/CLL_RIC/', pattern = '*.mgatk')) {
  boo = readRDS(paste0('./data/CLL_RIC/',library, '/final/mgatk.rds'))
  if (class(mgatk.data.RNA) == 'numeric') {
    mgatk.data.RNA = boo
  } else {
    mgatk.data.RNA = cbind(mgatk.data.RNA, boo)
  }
}
colnames(mgatk.data.RNA) = gsub(pattern = 'Tumor-', replacement = '', colnames(mgatk.data.RNA))
colnames(mgatk.data.RNA) = gsub(pattern = 'CLL-CRC-', replacement = '', colnames(mgatk.data.RNA))

# provide lower limit of detection for variant outside of indicated patient
detection_limit = function(variant, patient) {
  pos = as.numeric(gsub("\\D+", "", variant))
  REF = stringr::str_split_fixed(gsub("\\d", "", variant), pattern = ">", n = 2)[, 1]
  ALT = stringr::str_split_fixed(gsub("\\d", "", variant), pattern = ">", n = 2)[, 2]
  
  samples = meta$Sample[which(meta$Patient != patient)]
  
  refs <- paste0(REF, c("_counts_fw", "_counts_rev"))
  alts = paste0(ALT, c("_counts_fw", "_counts_rev"))
  
  REF.reads = 
  assays(mgatk.data.RNA)[[refs[1]]][as.numeric(pos), samples] +
    assays(mgatk.data.RNA)[[refs[2]]][as.numeric(pos), samples]
  
  ALT.reads = 
    assays(mgatk.data.RNA)[[alts[1]]][as.numeric(pos), samples] +
    assays(mgatk.data.RNA)[[alts[2]]][as.numeric(pos), samples]
  
  sum(ALT.reads) / sum(REF.reads)
}

### for each case plot mtDNA mutations with most extreme changes 
for (p in unique(meta$Patient)) {
  
  mtDNA.df = read.table(paste0('./data/CLL_RIC/processed_together/', p, '_de_novo.csv' ), 
                        check.names = F) %>%
    as.data.frame()
  rownames(mtDNA.df) = mtDNA.df$variant
  
  variants = mtDNA.df$variant[which(mtDNA.df$significant == 'significant')]
  # exclude variants that are not detectable at either end
  variants = variants[which(mtDNA.df[variants, "change"] != 0 & mtDNA.df[variants, "change"] != Inf)]
  
  # variants with change >10 or <0.1
  variants = variants[which(mtDNA.df[variants, "change"] > 2 | mtDNA.df[variants, "change"] < 0.5)]
  
  detection.limit.df = data.frame(variant = variants, 
                                  detection.limit = sapply(variants, FUN = function(x) { detection_limit(x, p)}))
  
  # filter variants with detection limit >0.1%
  variants = variants[which(detection.limit.df[variants, "detection.limit"] < 0.002)]
  if (isEmpty(variants)) {
    next()
  }
  
  detection.limit.df = detection.limit.df[which(detection.limit.df$detection.limit < 0.002),]
  
  timepoints = colnames(mtDNA.df)[grepl('heteroplasmy', colnames(mtDNA.df))]
  mtDNA.plot.data = mtDNA.df %>% dplyr::select(c('variant', timepoints)) %>%
    dplyr::filter(variant %in% variants) %>%
    pivot_longer(cols = timepoints)
  mtDNA.plot.data$day = sapply(mtDNA.plot.data$name, function(x) {
    meta$Day2[which(paste0(meta$Timepoint, '.heteroplasmy') == x & meta$Patient == p)]
  })
  mtDNA.plot.data$value.log = mtDNA.plot.data$value
  mtDNA.plot.data$value.log[which(mtDNA.plot.data$value < 0.0001)] = 0.0001
  
  color.scale = as.character(BuenColors::jdb_palette(name = 'corona', n = length(variants)))
  names(color.scale) = variants
  color.scale = c(color.scale, c('pre-HSCT' = 'blue', 'post-HSCT' = 'firebrick'))
  
  mtDNA.plot=ggplot() + 
    geom_hline(yintercept = 100*as.numeric(detection.limit.df[variants, "detection.limit"]), 
               color = color.scale[variants], linetype = 'dashed') +
    geom_line(data=mtDNA.plot.data, aes(x=day, y=100*value.log, 
                                        group=variant,color=as.character(variant))) +
    geom_point(data=mtDNA.plot.data, aes(x=day, y=100*value.log, 
                                         group=variant,color=as.character(variant))) +
    geom_point(data = meta[which(meta$Patient == p),], aes(x=Day2, y=250, color=Description), size=2) + 
    scale_y_log10('% heteroplasmy', limits = c(0.01,250), breaks = c(0.01, 0.1, 1, 10, 100), 
                  labels = c('0.01', '0.1', '1', '10', '100')) +
    scale_color_manual(values = color.scale) +
    ggtitle(p) + 
    theme_classic() +
    theme(legend.position = 'none',
          axis.title = element_text('Arial', size=10, color='black'),
          axis.text = element_text('Arial', size=10, color='black'),
          axis.title.x = element_blank(),
          plot.title = element_text('Arial', size=10, color='black', hjust = 0.5))
  ggsave(paste0('./figure_CLL_WES_RNA/figures/CLL_RIC/detection_limit/20240103_', p, '.svg'), 
         width = 1.5, height = 1.5, plot = mtDNA.plot)
}