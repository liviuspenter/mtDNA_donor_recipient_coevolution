# mixing across all 3 datasets
results.df.filtered.RNA = read.csv(file = './data/mixing/bulk/PRJNA741686.csv')
results.df.filtered.DNA = read.csv(file = './data/mixing/bulk/PRJNA486215_PRJNA563929.csv')
results.df.filtered = rbind(results.df.filtered.DNA, results.df.filtered.RNA)

# simulate data across all 3 datasets
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
      df = rbind(df, data.frame(individual1 = individual1, 
                                individual2 = individual2, 
                                variants = length(informative_variants),
                                dataset = ifelse(unique(results.df.filtered$dataset[which(results.df.filtered$sample == individual1)]) ==
                                                   unique(results.df.filtered$dataset[which(results.df.filtered$sample == individual2)]), 
                                                 'same', 'different')))
      
    }
  }
}

write.csv2(df, file = './data/mixing/bulk/mixing.all.csv', sep = '\t', quote = F)

ggplot(df, aes(x=variants, fill=dataset)) + geom_histogram(binwidth=1) +
  geom_vline(xintercept = 62) + 
  scale_fill_manual(values = c('different' = 'black', 'same' = 'grey')) + 
  scale_x_continuous('informative variants per pair', limits = c(-1,100)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/bulk/figures/20230327_informative_variants_all.svg', width = 2, height = 2)

# haplotype information from HaploGrep3
haplotypes = as.data.frame(data.table::fread('./data/mixing/bulk/all.haplotypes.csv'))
rownames(haplotypes) = paste0(haplotypes$SampleID, '.mgatk')
df$haplotype1 = haplotypes[df$individual1, 'Haplogroup']
df$haplotype2 = haplotypes[df$individual2, 'Haplogroup']

df$haplogroup1 = substr(df$haplotype1, start = 1, stop = 1)
df$haplogroup2 = substr(df$haplotype2, start = 1, stop = 1)

df$different.haplotype = factor(ifelse(df$haplogroup1 != df$haplogroup2, 'yes', 'no'), 
                                levels = c('yes', 'no'))

ggplot() + 
  geom_histogram(data=df[which(df$different.haplotype == 'yes'),], 
                 aes(x=variants, fill=different.haplotype, alpha=different.haplotype), binwidth=1) +
  geom_histogram(data=df[which(df$different.haplotype == 'no'),], 
                 aes(x=variants, fill=different.haplotype, alpha=different.haplotype), binwidth=1) +
  scale_fill_manual(values = c('yes' = 'grey', 'no' = 'orange')) + 
  scale_x_continuous('informative variants per pair', limits = c(-1,100)) +
  scale_alpha_manual(values = c('yes' = 1, 'no' = 0.7)) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/bulk/figures/20230711_informative_variants_all_haplogroup.svg', width = 2, height = 2)
