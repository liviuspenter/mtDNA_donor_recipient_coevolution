setwd('/Users/shaka87/dfci/asap_seq/')

# Gene locations
y.max = 20
GenePos.tib <- tibble(Names = c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3",
                                "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
  mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(1.3*y.max,
                                                            1.2*y.max,
                                                            1.1*y.max), length.out = 15))

GenePos.tib$Names = stringr::str_split_fixed(GenePos.tib$Names, pattern = '\\.', n=2)[,2]

germline.variants = read.csv2('./data/objects/20220117_IST_germline_variants.csv')
germline.variants = rbind(germline.variants, read.csv2('./data/objects/20220117_AML_germline_variants.csv'))

germline.variants$position = gsub("([0-9]+).*$", "\\1", germline.variants$variant)

ggplot(germline.variants, aes(x=as.numeric(position))) + geom_histogram(col='black', binwidth=100) +
  geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
  geom_text(data = GenePos.tib, aes(x = mid, y = 1.03*ycoord, label = Names), size = 2) +
  scale_x_continuous('chrM', breaks = seq(0,16000, 2000)) + 
  scale_y_continuous('# mtDNA variants') + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=8, color='black'),
        axis.text = element_text('Arial', size=8, color='black'))
ggsave('./figures/plots/20220119_mtDNA_germline_variants.svg', width = 3.5, height = 2.5)

boo = as.data.frame(table(germline.variants$position))

ggplot(boo, aes(x=as.numeric(unfactor(Var1)), y=Freq)) + geom_col() + theme_classic()