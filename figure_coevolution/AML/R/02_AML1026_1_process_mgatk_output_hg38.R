# process mgatk output to perform mtDNA variant calling and visualize mtDNA coverage

library(ggplot2)
library(Seurat)
library(Signac)

AML1026.1.mito = ReadMGATK('./data/coevolution/1026_1.mgatk/')
AML1026.1.so <- CreateAssayObject(counts = AML1026.1.mito$counts)
AML1026.1.mgatk = CreateSeuratObject(counts = AML1026.1.so, project = 'AML1026.1.mgatk', assay = 'mito')
rm(AML1026.1.so)
AML1026.1.mgatk = AddMetaData(AML1026.1.mgatk, metadata = AML1026.1.mito$depth, col.name = 'mtDNA_depth')

variable.sites <- IdentifyVariants(AML1026.1.mgatk, assay = "mito", refallele = AML1026.1.mito$refallele)

VariantPlot(variants = variable.sites)
high.conf = variable.sites[which(variable.sites$n_cells_conf_detected > 5),]

AML1026.1.mgatk <- AlleleFreq(object = AML1026.1.mgatk, variants = high.conf$variant, assay = "mito")
AML1026.1.vaf = as.data.frame(t(as.data.frame(GetAssayData(AML1026.1.mgatk[["alleles"]]))))
write.csv(AML1026.1.vaf, file = './data/coevolution/vafs/20220211_AML1026.1_mtDNA_vafs.csv')
write.table(variable.sites, file = './data/coevolution/vafs/20220211_AML1026_1_variable_sites.csv')

variable.sites = read.csv2(file = './data/coevolution/vafs/20220211_AML1026_1_variable_sites.csv', sep = ' ')

# plot coverage
GenePos.tib <- tibble(Names = c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3",
                                "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
  mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(280,300,320), length.out = 15))

p=ggplot(variable.sites, aes(x=as.numeric(position), y=as.numeric(mean_coverage))) + geom_col() +
  geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
  geom_text(data = GenePos.tib, aes(x = mid, y = ycoord+8, label = Names), size = 2) +
  scale_x_continuous('chrM') +
  scale_y_continuous('mean coverage') +
  theme_classic() +
  theme(axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./figure_coevolution/AML/figures/20220211_AML1026_1_coverage.png', width = 3.5, height = 2.5, plot = p, dpi = 600)
