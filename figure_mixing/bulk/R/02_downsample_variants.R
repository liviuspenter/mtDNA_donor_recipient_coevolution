# downsample experiment
results.df.filtered.RNA = read.csv(file = './data/mixing/bulk/PRJNA741686.csv')
results.df.filtered.DNA = read.csv(file = './data/mixing/bulk/PRJNA486215_PRJNA563929.csv')
results.df.filtered = rbind(results.df.filtered.DNA, results.df.filtered.RNA)

samples = unique(results.df.filtered$sample)
df = data.frame()
for (i in seq(1,length(samples))) {
  boo = results.df.filtered %>% filter(sample %in% samples[seq(1,i)])
  number_of_variants = length(unique(boo$variant))
  df = rbind(df, data.frame(n = i, variants = number_of_variants))
}
ggplot(df, aes(x=n, y=variants)) + 
  geom_line(color = 'grey') +
  geom_point(size=0.5) + 
  scale_x_continuous('individuals') + 
  scale_y_continuous('mtDNA variants') + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/bulk/figures/20230319_mtDNA_variants_downsampling.svg', width = 2, height = 2)