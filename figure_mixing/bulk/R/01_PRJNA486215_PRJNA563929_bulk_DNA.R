setwd('/Users/shaka87/dfci/asap_seq/')

library(dplyr)
library(ggplot2)

source('./analysis/bulk_mtDNA/00_variant_calling_atac.R')

GenePos.tib <- tibble(Names = c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3",
                                "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))

### plot coverage
coverage.df = data.frame()
for (library in list.files('./data/bulk/PRJNA486215/', pattern = '*.mgatk')) {
  message(library)
  coverage = data.table::fread(paste0('./data/bulk/PRJNA486215/',library, '/final/mgatk.coverage.txt.gz'))
  colnames(coverage) = c('pos','sample', 'coverage')
  coverage.df = rbind(coverage.df, coverage)
}
coverage = coverage.df %>% group_by(pos) %>% summarize(coverage.mean = mean(coverage))
for (library in list.files('./data/bulk/PRJNA563929/', pattern = '*.mgatk')) {
  message(library)
  coverage = data.table::fread(paste0('./data/bulk/PRJNA563929/',library, '/final/mgatk.coverage.txt.gz'))
  colnames(coverage) = c('pos','sample', 'coverage')
  coverage.df = rbind(coverage.df, coverage)
}
coverage = coverage.df %>% group_by(pos) %>% summarize(coverage.mean = mean(coverage))

y.max = max(coverage$coverage.mean)/1000
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
  mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(1.3*y.max,
                                                            1.2*y.max,
                                                            1.1*y.max), length.out = 15))

for (i in seq(1, nrow(GenePos.tib))) {
  coverage[which(between(coverage$pos, GenePos.tib$start[i], GenePos.tib$end[i])), 'gene'] = GenePos.tib$Names[i]
}

ggplot(coverage, aes(x=pos, y=coverage.mean/1000)) + 
  ggrastr::rasterize(geom_col(aes(fill=gene)), dpi=600) +
  geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord, color=Names)) +
  geom_text(data = GenePos.tib, aes(x = mid, y = 1.03*ycoord, label = Names), size = 2) +
  scale_x_continuous('chrM', breaks = seq(0,16000, 4000)) + 
  scale_y_continuous('1000x coverage') + 
  scale_fill_manual(values = BuenColors::jdb_palette(name = 'corona', n=15), na.value = 'black') +
  scale_color_manual(values = BuenColors::jdb_palette(name = 'corona', n=15), na.value = 'black') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figures/bulk/20230319_mtDNA_coverage_bulk_DNA.svg', width = 3, height = 2.5)



### read variants
results.df = data.frame()
for (library in list.files('./data/bulk/PRJNA486215/', pattern = '*.mgatk')) {
  message(library)
  results = variant_calling_bulk(sample.bulk = readRDS(paste0('./data/bulk/PRJNA486215/',library, '/final/mgatk.rds')), 
                                 coverage.position = 20, strand.coordination = 0, total.coverage = 30,
                                 sample.name = s)
  results = results[which(!is.na(results$heteroplasmy)),]
  results$sample = library
  results$dataset = 'PRJNA486215'
  
  results.df = rbind(results.df, results)
}
for (library in list.files('./data/bulk/PRJNA563929//', pattern = '*.mgatk')) {
  message(library)
  results = variant_calling_bulk(sample.bulk = readRDS(paste0('./data/bulk/PRJNA563929/',library, '/final/mgatk.rds')), 
                                 coverage.position = 20, strand.coordination = 0, total.coverage = 30,
                                 sample.name = s)
  results = results[which(!is.na(results$heteroplasmy)),]
  results$sample = library
  results$dataset = 'PRJNA563929'
  
  results.df = rbind(results.df, results)
}
results.df$pos = as.numeric(gsub('\\D+','', results.df$variant))
for (i in seq(1, nrow(GenePos.tib))) {
  results.df[which(between(results.df$pos, GenePos.tib$start[i], GenePos.tib$end[i])), 'gene'] = GenePos.tib$Names[i]
}

# only consider variants with heteroplasmy >98%
results.df.filtered = results.df %>% filter(heteroplasmy > 0.98)

### variants per patient
write.csv(results.df.filtered, file = './data/bulk/PRJNA486215_PRJNA563929.csv', quote = F, row.names = F)
results.df.filtered = read.csv(file = './data/bulk/PRJNA486215_PRJNA563929.csv')

### visualize frequency across chrM
df = as.data.frame(table(results.df.filtered$variant))
df$percentage = df$Freq / length(unique(results.df$sample))
df$pos = as.numeric(gsub('\\D+','', df$Var1))
df$variant = gsub('[[:digit:]]', '', df$Var1)

# Gene locations
y.max = 30
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
  mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(1.3*y.max,
                                                            1.2*y.max,
                                                            1.1*y.max), length.out = 15))

df = df %>% mutate(pos.bin = cut(pos, breaks = seq(1,16600,100), labels=F))
boo = df %>% group_by(pos.bin) %>% tally()
boo$pos = 100*boo$pos.bin
for (i in seq(1, nrow(GenePos.tib))) {
  boo[which(between(boo$pos, GenePos.tib$start[i], GenePos.tib$end[i])), 'gene'] = GenePos.tib$Names[i]
}

ggplot(boo, aes(x=pos, y=n)) + geom_col(aes(fill=gene)) +
  geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord, color=Names)) +
  geom_text(data = GenePos.tib, aes(x = mid, y = 1.03*ycoord, label = Names), size = 2) +
  scale_x_continuous('chrM', breaks = seq(0,16000, 4000), limits = c(0,17000)) + 
  scale_y_continuous('# mtDNA variants') + 
  scale_color_manual(values = BuenColors::jdb_palette(name = 'corona', n=15), na.value = 'black') +
  scale_fill_manual(values = BuenColors::jdb_palette(name = 'corona', n=15), na.value = 'black') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figures/bulk/20230319_mtDNA_variants_bulk_DNA.svg', width = 3, height = 2.5)

### calculate number of variants for each pair

# plot variants per sample
df = results.df.filtered %>% group_by(sample) %>% summarize(n = length(variant))
# median 30 variants
median(df$n)
# max 78 variants
max(df$n)
# min 7 variants 
min(df$n)

ggplot(df, aes(x=n)) + geom_histogram(binwidth = 1, fill='firebrick') +
  scale_x_continuous('mtDNA variants per individual', limits = c(-1,100)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figures/bulk/20230319_absolute_mtDNA_variants_bulk_DNA.svg', width = 2, height = 2)

# simulate data
df = data.frame()
informative_variants_all = c()
for (individual1 in unique(results.df.filtered$sample)) {
  message(individual1)
  for (individual2 in unique(results.df.filtered$sample)) {
    if (individual1 != individual2) {
      informative_variants = c(setdiff(results.df.filtered$variant[which(results.df.filtered$sample == individual1)],
                                     results.df.filtered$variant[which(results.df.filtered$sample == individual2)]),
                               setdiff(results.df.filtered$variant[which(results.df.filtered$sample == individual2)],
                                       results.df.filtered$variant[which(results.df.filtered$sample == individual1)]))
      informative_variants_all = c(informative_variants_all, informative_variants)
      df = rbind(df, data.frame(individual1 = individual1, individual2 = individual2, variants = length(informative_variants)))
      
    }
  }
}

# plot number of informative variants
ggplot(df, aes(x=variants)) + geom_histogram(binwidth = 1, fill='firebrick') +
  geom_vline(xintercept = 62) + 
  scale_x_continuous('informative variants per pair', limits = c(-1,100)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figures/bulk/20230319_informative_mtDNA_variants_bulk_DNA.svg', width = 2, height = 2)

# median 35
median(df$variants)
# max 94
max(df$variants)
# 4.9% high variance
length(which(df$variants > 62)) / nrow(df)
# 0% no variant 
length(which(df$variants == 0)) / nrow(df)

# plot location of informative variants
df =  as.data.frame(as.numeric(gsub('\\D+','', informative_variants_all)) )
colnames(df) = 'pos'
df = df %>% mutate(pos.bin = cut(pos, breaks = seq(1,16600,100), labels=F))
boo = df %>% group_by(pos.bin) %>% tally()
boo$pos = 100*boo$pos.bin
for (i in seq(1, nrow(GenePos.tib))) {
  boo[which(between(boo$pos, GenePos.tib$start[i], GenePos.tib$end[i])), 'gene'] = GenePos.tib$Names[i]
}

ggplot(data=boo, aes(x=pos, y=n)) + geom_col(aes(fill=gene)) + 
  scale_x_continuous('chrM', breaks = seq(0,16000, 4000)) + 
  scale_y_continuous('times used for deconvolution') + 
  scale_fill_manual(values = BuenColors::jdb_palette(name = 'corona', n=15), na.value = 'black') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figures/bulk/20230319_informative_mtDNA_variants_distribution_bulk_DNA.svg', width = 3, height = 2)
