SNP.count = data.table::fread('./data/mixing/bulk/PRJNA741686/SNP_count.csv')
ggplot(SNP.count, aes(x=V2)) + geom_histogram(binwidth = 10000, fill='lightblue') + 
  scale_x_continuous('expressed SNPs per individual', limits = c(0,150000), breaks = c(0,50000,100000), 
                     labels = c("0","50000","100000")) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/bulk/figures/20230330_expressed_SNPs_per_individual.svg', width = 2, height = 2)

SNP.data = data.table::fread('./data/mixing/bulk/PRJNA741686/simulated_snps.csv.gz')
SNP.data$sum = SNP.data$V3 + SNP.data$V4

ggplot(SNP.data, aes(x=sum)) + geom_histogram(binwidth = 500, fill='lightblue') + 
  scale_x_continuous('informative SNPs per pair', limits = c(0, 150000), breaks = c(0,50000,100000), 
                     labels = c("0","50000","100000")) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_mixing/bulk/figures/20230330_informative_SNPs_per_pair.svg', width = 2, height = 2)