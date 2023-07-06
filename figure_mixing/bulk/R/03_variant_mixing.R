# downsample experiment
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

ggplot(df, aes(x=variants, fill=dataset)) + geom_histogram(binwidth=1) +
  geom_vline(xintercept = 62) + 
  scale_fill_manual(values = c('different' = 'black', 'same' = 'grey')) + 
  scale_x_continuous('informative variants per pair', limits = c(-1,100)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/bulk/figures/20230327_informative_variants_all.svg', width = 2, height = 2)
